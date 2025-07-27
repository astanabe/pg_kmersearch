#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use Time::HiRes qw(time);

# PostgreSQL connection parameters
my $dbname = $ENV{PGDATABASE} || 'postgres';
my $host = $ENV{PGHOST} || 'localhost';
my $port = $ENV{PGPORT} || '5432';
my $user = $ENV{PGUSER} || $ENV{USER};

# Connect to database
my $dbh = DBI->connect("dbi:Pg:dbname=$dbname;host=$host;port=$port", $user, '', 
    {AutoCommit => 1, RaiseError => 1})
    or die "Cannot connect to database: $DBI::errstr\n";

# SIMD capabilities for AMD64
my @simd_capabilities = (
    { name => 'AVX512VBMI2', value => 6 },
    { name => 'AVX512VBMI',  value => 5 },
    { name => 'AVX512BW',    value => 4 },
    { name => 'AVX512F',     value => 3 },
    { name => 'BMI2',        value => 2 },
    { name => 'AVX2',        value => 1 },
    { name => 'None',        value => 0 },
);

# Test DNA sequences
my @test_sequences = (
    # Short sequences
    'ACGTACGTACGT',
    'TGCATGCATGCA',
    # Medium sequences
    'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT',
    'TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA',
    # Long sequences
    'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT' x 10,
    'TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA' x 10,
    # Sequences with all bases
    'AAAACCCCGGGGTTTT' x 8,
    'ACGTACGTACGTACGT' x 8,
    # Random-like sequences
    'ACGTAGCTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT',
    'TGCATCGATGCATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA',
);

# Test DNA4 sequences with degenerate bases
my @test_dna4_sequences = (
    'ACGTNNNNACGT',
    'MRWSYKVHDBN',
    'ACGTMRWSACGT',
    'NNNNNNNNNNNN',
    'ACGTACGTACGTMRWSYKVHDBNACGTACGTACGT',
);

# Create test table
print "Creating test table...\n";
$dbh->do("DROP TABLE IF EXISTS simd_test CASCADE");
$dbh->do("CREATE TABLE simd_test (id serial PRIMARY KEY, dna2_seq dna2, dna4_seq dna4)");

# Prepare insert statement
my $insert_sth = $dbh->prepare("INSERT INTO simd_test (dna2_seq, dna4_seq) VALUES (?, ?)");

# Get initial SIMD capability
my ($auto_capability) = $dbh->selectrow_array("SELECT kmersearch_simd_capability()");
print "Auto-detected SIMD capability: $auto_capability\n\n";

# Test each SIMD level
foreach my $simd (@simd_capabilities) {
    print "=" x 60 . "\n";
    print "Testing with SIMD capability: $simd->{name} (value: $simd->{value})\n";
    print "=" x 60 . "\n";
    
    # Set SIMD capability
    $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
    
    # Verify setting
    my ($current) = $dbh->selectrow_array("SELECT kmersearch_simd_capability()");
    print "Current SIMD capability: $current\n";
    
    # Clear table
    $dbh->do("TRUNCATE TABLE simd_test");
    
    # Test DNA2 encoding/decoding
    print "\nTesting DNA2 encoding/decoding...\n";
    my $dna2_start = time();
    my $dna2_count = 0;
    
    foreach my $seq (@test_sequences) {
        # Insert
        $insert_sth->execute($seq, undef);
        $dna2_count++;
        
        # Read back and verify
        my ($retrieved) = $dbh->selectrow_array(
            "SELECT dna2_seq::text FROM simd_test WHERE id = ?", 
            undef, $dbh->last_insert_id(undef, undef, 'simd_test', 'id')
        );
        
        if ($retrieved ne $seq) {
            die "DNA2 mismatch! Original: $seq, Retrieved: $retrieved\n";
        }
    }
    
    # Batch operations
    for (my $i = 0; $i < 100; $i++) {
        foreach my $seq (@test_sequences) {
            $insert_sth->execute($seq, undef);
            $dna2_count++;
        }
    }
    
    my $dna2_time = time() - $dna2_start;
    printf "  Processed %d DNA2 sequences in %.3f seconds (%.0f seq/sec)\n", 
           $dna2_count, $dna2_time, $dna2_count / $dna2_time;
    
    # Test DNA4 encoding/decoding
    print "\nTesting DNA4 encoding/decoding...\n";
    my $dna4_start = time();
    my $dna4_count = 0;
    
    foreach my $seq (@test_dna4_sequences) {
        # Insert
        $insert_sth->execute(undef, $seq);
        $dna4_count++;
        
        # Read back and verify
        my ($retrieved) = $dbh->selectrow_array(
            "SELECT dna4_seq::text FROM simd_test WHERE id = ?", 
            undef, $dbh->last_insert_id(undef, undef, 'simd_test', 'id')
        );
        
        if ($retrieved ne $seq) {
            die "DNA4 mismatch! Original: $seq, Retrieved: $retrieved\n";
        }
    }
    
    # Batch operations
    for (my $i = 0; $i < 100; $i++) {
        foreach my $seq (@test_dna4_sequences) {
            $insert_sth->execute(undef, $seq);
            $dna4_count++;
        }
    }
    
    my $dna4_time = time() - $dna4_start;
    printf "  Processed %d DNA4 sequences in %.3f seconds (%.0f seq/sec)\n", 
           $dna4_count, $dna4_time, $dna4_count / $dna4_time;
    
    # Test k-mer operations
    print "\nTesting k-mer extraction...\n";
    $dbh->do("SET kmersearch.kmer_size = 8");
    
    my $kmer_start = time();
    my $kmer_count = 0;
    
    # Create GIN index
    $dbh->do("CREATE INDEX simd_test_dna2_gin ON simd_test USING gin(dna2_seq gin_kmer_ops)");
    
    # Perform searches
    foreach my $query (@test_sequences[0..4]) {
        my $count = $dbh->selectrow_array(
            "SELECT COUNT(*) FROM simd_test WHERE dna2_seq =% ?::dna2",
            undef, $query
        );
        $kmer_count += $count;
    }
    
    my $kmer_time = time() - $kmer_start;
    printf "  Performed k-mer searches in %.3f seconds\n", $kmer_time;
    
    # Drop index for next iteration
    $dbh->do("DROP INDEX simd_test_dna2_gin");
    
    print "\n";
}

# Reset to auto-detect
print "Resetting to auto-detected SIMD capability...\n";
$dbh->do("SET kmersearch.force_simd_capability = -1");
my ($final) = $dbh->selectrow_array("SELECT kmersearch_simd_capability()");
print "Final SIMD capability: $final\n";

# Cleanup
$dbh->do("DROP TABLE simd_test CASCADE");
$dbh->disconnect();

print "\nAll tests completed successfully!\n";