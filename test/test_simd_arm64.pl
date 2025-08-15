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
    {AutoCommit => 1, RaiseError => 1, PrintError => 1})
    or die "Cannot connect to database: $DBI::errstr\n";

# Create extension if not exists
print "Creating/verifying pg_kmersearch extension...\n";
eval {
    $dbh->do("CREATE EXTENSION IF NOT EXISTS pg_kmersearch");
};
if ($@) {
    die "Failed to create pg_kmersearch extension: $@\n";
}
print "Extension pg_kmersearch is ready.\n\n";

# Test parameters
my $sequence_length = 4096;
my $total_sequences = 10000;
my $num_query_sequences = 1000;
my $max_appearance_rate = 0.1;
my $kmer_size = 16;

# SIMD capabilities for ARM64 - using numeric values for force_simd_capability
my @simd_capabilities = (
    { name => 'None',    value => 0 },
    { name => 'NEON',    value => 21 },
    { name => 'SVE',     value => 22 },
    { name => 'SVE2',    value => 23 },
);

# Generate random DNA sequences (shared across all tests)
print "Generating $total_sequences random DNA sequences of $sequence_length bp...\n";
my @sequences;
my @bases = ('A', 'C', 'G', 'T');

for (my $i = 0; $i < $total_sequences; $i++) {
    my $seq = '';
    for (my $j = 0; $j < $sequence_length; $j++) {
        $seq .= $bases[int(rand(4))];
    }
    push @sequences, $seq;
    if (($i + 1) % 1000 == 0) {
        print "  Generated " . ($i + 1) . " sequences...\n";
    }
}

# Generate DNA4 sequences with degenerate codes
print "Generating $total_sequences random DNA4 sequences with degenerate codes...\n";
my @dna4_sequences;
my @dna4_bases = ('A', 'C', 'G', 'T', 'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B', 'N');

for (my $i = 0; $i < $total_sequences; $i++) {
    my $seq = '';
    for (my $j = 0; $j < $sequence_length; $j++) {
        # Use mostly standard bases with very occasional degenerate codes
        if (rand() < 0.001) {  # 0.1% chance for degenerate code
            $seq .= $dna4_bases[int(rand(@dna4_bases))];
        } else {
            $seq .= $bases[int(rand(4))];
        }
    }
    push @dna4_sequences, $seq;
    if (($i + 1) % 1000 == 0) {
        print "  Generated " . ($i + 1) . " DNA4 sequences...\n";
    }
}

# Select query sequences (first 1000)
my @query_sequences = @sequences[0..($num_query_sequences-1)];
my @dna4_query_sequences = @dna4_sequences[0..($num_query_sequences-1)];

print "Generated " . scalar(@sequences) . " DNA2 sequences\n";
print "Generated " . scalar(@dna4_sequences) . " DNA4 sequences\n";
print "Selected " . scalar(@query_sequences) . " DNA2 query sequences\n";
print "Selected " . scalar(@dna4_query_sequences) . " DNA4 query sequences\n\n";

# Get initial SIMD capability
my ($auto_capability) = $dbh->selectrow_array("SELECT kmersearch_simd_capability()");
print "Auto-detected SIMD capability: $auto_capability\n\n";

# Store results for comparison
my %test_results;
my %timing_results;
my %additional_io_timing;

# Function to measure execution time with error handling
sub measure_time {
    my ($name, $code) = @_;
    my $start = time();
    my $result;
    eval {
        local $SIG{ALRM} = sub { die "Timeout in $name\n" };
        alarm(300);  # 5 minute timeout
        $result = $code->();
        alarm(0);
    };
    my $end = time();
    if ($@) {
        print "ERROR in $name: $@\n";
        return ($end - $start, undef);
    }
    return ($end - $start, $result);
}

# Create base table function
sub create_base_table {
    my ($table_name, $data_type) = @_;
    $data_type ||= 'DNA2';  # Default to DNA2
    
    # Drop existing table
    $dbh->do("DROP TABLE IF EXISTS $table_name CASCADE");
    
    # Create table
    $dbh->do("CREATE TABLE $table_name (id SERIAL PRIMARY KEY, seq $data_type)");
    
    # Insert sequences
    my $sth = $dbh->prepare("INSERT INTO $table_name (seq) VALUES (?)");
    my $sequences_ref = ($data_type eq 'DNA4') ? \@dna4_sequences : \@sequences;
    foreach my $seq (@$sequences_ref) {
        $sth->execute($seq);
    }
    $sth->finish();
}

# Test execution order: DNA2 scalar I/O -> DNA4 scalar I/O -> DNA2 scalar input/SIMD output -> ...
print "=" x 70 . "\n";
print "SIMD Implementation Verification Tests\n";
print "Following test order: DNA2 scalar -> DNA4 scalar -> ...\n";
print "=" x 70 . "\n\n";

# DNA2 Test 1: Scalar input/Scalar output
print "DNA2 Test 1: Scalar input (force_simd_capability=0) / Scalar output (force_simd_capability=0)\n";
my ($dna2_scalar_io_time, $dna2_scalar_io_result) = measure_time("dna2_scalar_io", sub {
    my $table_name = "simd_test_dna2_scalar_io";
    
    # Set to scalar mode for input
    $dbh->do("SET kmersearch.force_simd_capability = 0");
    
    # Create table and insert data
    create_base_table($table_name);
    
    # Verify with scalar output
    $dbh->do("SET kmersearch.force_simd_capability = 0");
    
    my $sth = $dbh->prepare("SELECT id, seq::text FROM $table_name ORDER BY id");
    $sth->execute();
    
    my $match_count = 0;
    my $row_count = 0;
    while (my ($id, $seq) = $sth->fetchrow_array()) {
        if ($sequences[$id-1] eq $seq) {
            $match_count++;
        } else {
            print "  Mismatch at ID $id: expected '$sequences[$id-1]', got '$seq'\n" if $row_count < 5;
        }
        $row_count++;
    }
    $sth->finish();
    
    # Clean up
    $dbh->do("DROP TABLE $table_name");
    
    return {
        total_rows => $row_count,
        matched => $match_count,
        success => ($match_count == scalar(@sequences))
    };
});

$additional_io_timing{dna2_scalar_io} = $dna2_scalar_io_time;
printf "  Time: %.3fs, Result: %d/%d sequences match (%.1f%%)\n", 
       $dna2_scalar_io_time, $dna2_scalar_io_result->{matched}, $dna2_scalar_io_result->{total_rows},
       ($dna2_scalar_io_result->{matched} / $dna2_scalar_io_result->{total_rows}) * 100;

# DNA4 Test 1: Scalar input/Scalar output
print "\nDNA4 Test 1: Scalar input (force_simd_capability=0) / Scalar output (force_simd_capability=0)\n";
my ($dna4_scalar_io_time, $dna4_scalar_io_result) = measure_time("dna4_scalar_io", sub {
    my $table_name = "simd_test_dna4_scalar_io";
    
    # Set to scalar mode for input
    $dbh->do("SET kmersearch.force_simd_capability = 0");
    
    # Create table and insert data
    create_base_table($table_name, 'DNA4');
    
    # Verify with scalar output
    $dbh->do("SET kmersearch.force_simd_capability = 0");
    
    my $sth = $dbh->prepare("SELECT id, seq::text FROM $table_name ORDER BY id");
    $sth->execute();
    
    my $match_count = 0;
    my $row_count = 0;
    while (my ($id, $seq) = $sth->fetchrow_array()) {
        if ($dna4_sequences[$id-1] eq $seq) {
            $match_count++;
        } else {
            print "  Mismatch at ID $id: expected '$dna4_sequences[$id-1]', got '$seq'\n" if $row_count < 5;
        }
        $row_count++;
    }
    $sth->finish();
    
    # Clean up
    $dbh->do("DROP TABLE $table_name");
    
    return {
        total_rows => $row_count,
        matched => $match_count,
        success => ($match_count == scalar(@dna4_sequences))
    };
});

$additional_io_timing{dna4_scalar_io} = $dna4_scalar_io_time;
printf "  Time: %.3fs, Result: %d/%d sequences match (%.1f%%)\n", 
       $dna4_scalar_io_time, $dna4_scalar_io_result->{matched}, $dna4_scalar_io_result->{total_rows},
       ($dna4_scalar_io_result->{matched} / $dna4_scalar_io_result->{total_rows}) * 100;

# DNA2 Test 2: Scalar input/SIMD output (all SIMD values)
print "\nDNA2 Test 2: Scalar input (force_simd_capability=0) / SIMD output (all force_simd_capability>0)\n";
my ($dna2_scalar_input_time, $dna2_scalar_input_result) = measure_time("dna2_scalar_input_simd_output", sub {
    my $table_name = "simd_test_scalar_input";
    my %results;
    
    # Set to scalar mode for input
    $dbh->do("SET kmersearch.force_simd_capability = 0");
    
    # Create table and insert data
    create_base_table($table_name);
    
    # Test with each SIMD capability for output
    foreach my $simd (@simd_capabilities) {
        next if $simd->{value} == 0;  # Skip scalar
        
        my $test_start = time();
        $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
        
        my $sth = $dbh->prepare("SELECT id, seq::text FROM $table_name ORDER BY id");
        $sth->execute();
        
        my $match_count = 0;
        my $row_count = 0;
        while (my ($id, $seq) = $sth->fetchrow_array()) {
            if ($sequences[$id-1] eq $seq) {
                $match_count++;
            } else {
                print "  Mismatch with $simd->{name} output at ID $id\n" if $row_count < 5;
            }
            $row_count++;
        }
        $sth->finish();
        my $test_time = time() - $test_start;
        
        $results{$simd->{name}} = {
            time => $test_time,
            matched => $match_count,
            total => $row_count,
            success_rate => ($match_count / $row_count) * 100
        };
        
        printf "  $simd->{name} output: %d/%d sequences match (%.1f%%) - Time: %.3fs\n", 
               $match_count, $row_count, ($match_count / $row_count) * 100, $test_time;
    }
    
    # Clean up
    $dbh->do("DROP TABLE $table_name");
    
    return \%results;
});

$additional_io_timing{dna2_scalar_input_simd_output} = $dna2_scalar_input_time;

# DNA4 Test 2: Scalar input/SIMD output (all SIMD values)
print "\nDNA4 Test 2: Scalar input (force_simd_capability=0) / SIMD output (all force_simd_capability>0)\n";
my ($dna4_scalar_input_time, $dna4_scalar_input_result) = measure_time("dna4_scalar_input_simd_output", sub {
    my $table_name = "simd_test_dna4_scalar_input";
    my %results;
    
    # Set to scalar mode for input
    $dbh->do("SET kmersearch.force_simd_capability = 0");
    
    # Create table and insert data
    create_base_table($table_name, 'DNA4');
    
    # Test with each SIMD capability for output
    foreach my $simd (@simd_capabilities) {
        next if $simd->{value} == 0;  # Skip scalar
        
        my $test_start = time();
        $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
        
        my $sth = $dbh->prepare("SELECT id, seq::text FROM $table_name ORDER BY id");
        $sth->execute();
        
        my $match_count = 0;
        my $row_count = 0;
        while (my ($id, $seq) = $sth->fetchrow_array()) {
            if ($dna4_sequences[$id-1] eq $seq) {
                $match_count++;
            } else {
                print "  Mismatch with $simd->{name} output at ID $id\n" if $row_count < 5;
            }
            $row_count++;
        }
        $sth->finish();
        my $test_time = time() - $test_start;
        
        $results{$simd->{name}} = {
            time => $test_time,
            matched => $match_count,
            total => $row_count,
            success_rate => ($match_count / $row_count) * 100
        };
        
        printf "  $simd->{name} output: %d/%d sequences match (%.1f%%) - Time: %.3fs\n", 
               $match_count, $row_count, ($match_count / $row_count) * 100, $test_time;
    }
    
    # Clean up
    $dbh->do("DROP TABLE $table_name");
    
    return \%results;
});

$additional_io_timing{dna4_scalar_input_simd_output} = $dna4_scalar_input_time;

# DNA2 Test 3: SIMD input/Scalar output (all SIMD values)
print "\nDNA2 Test 3: SIMD input (all force_simd_capability>0) / Scalar output (force_simd_capability=0)\n";
my ($dna2_simd_input_time, $dna2_simd_input_result) = measure_time("dna2_simd_input_scalar_output", sub {
    my %results;
    
    foreach my $simd (@simd_capabilities) {
        next if $simd->{value} == 0;  # Skip scalar
        
        my $test_start = time();
        my $table_name = "simd_test_" . $simd->{name} . "_input";
        $table_name =~ s/[^a-zA-Z0-9_]/_/g;  # Replace non-alphanumeric chars
        
        # Set SIMD mode for input
        $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
        
        # Create table and insert data
        create_base_table($table_name);
        
        # Set to scalar mode for output
        $dbh->do("SET kmersearch.force_simd_capability = 0");
        
        my $sth = $dbh->prepare("SELECT id, seq::text FROM $table_name ORDER BY id");
        $sth->execute();
        
        my $match_count = 0;
        my $row_count = 0;
        while (my ($id, $seq) = $sth->fetchrow_array()) {
            if ($sequences[$id-1] eq $seq) {
                $match_count++;
            } else {
                print "  Mismatch with $simd->{name} input at ID $id\n" if $row_count < 5;
            }
            $row_count++;
        }
        $sth->finish();
        
        # Clean up
        $dbh->do("DROP TABLE $table_name");
        
        my $test_time = time() - $test_start;
        
        $results{$simd->{name}} = {
            time => $test_time,
            matched => $match_count,
            total => $row_count,
            success_rate => ($match_count / $row_count) * 100
        };
        
        printf "  $simd->{name} input: %d/%d sequences match (%.1f%%) - Time: %.3fs\n", 
               $match_count, $row_count, ($match_count / $row_count) * 100, $test_time;
    }
    
    return \%results;
});

$additional_io_timing{dna2_simd_input_scalar_output} = $dna2_simd_input_time;

# DNA4 Test 3: SIMD input/Scalar output (all SIMD values)
print "\nDNA4 Test 3: SIMD input (all force_simd_capability>0) / Scalar output (force_simd_capability=0)\n";
my ($dna4_simd_input_time, $dna4_simd_input_result) = measure_time("dna4_simd_input_scalar_output", sub {
    my %results;
    
    foreach my $simd (@simd_capabilities) {
        next if $simd->{value} == 0;  # Skip scalar
        
        my $test_start = time();
        my $table_name = "simd_test_dna4_" . $simd->{name} . "_input";
        $table_name =~ s/[^a-zA-Z0-9_]/_/g;  # Replace non-alphanumeric chars
        
        # Set SIMD mode for input
        $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
        
        # Create table and insert data
        create_base_table($table_name, 'DNA4');
        
        # Set to scalar mode for output
        $dbh->do("SET kmersearch.force_simd_capability = 0");
        
        my $sth = $dbh->prepare("SELECT id, seq::text FROM $table_name ORDER BY id");
        $sth->execute();
        
        my $match_count = 0;
        my $row_count = 0;
        while (my ($id, $seq) = $sth->fetchrow_array()) {
            if ($dna4_sequences[$id-1] eq $seq) {
                $match_count++;
            } else {
                print "  Mismatch with $simd->{name} input at ID $id\n" if $row_count < 5;
            }
            $row_count++;
        }
        $sth->finish();
        
        # Clean up
        $dbh->do("DROP TABLE $table_name");
        
        my $test_time = time() - $test_start;
        
        $results{$simd->{name}} = {
            time => $test_time,
            matched => $match_count,
            total => $row_count,
            success_rate => ($match_count / $row_count) * 100
        };
        
        printf "  $simd->{name} input: %d/%d sequences match (%.1f%%) - Time: %.3fs\n", 
               $match_count, $row_count, ($match_count / $row_count) * 100, $test_time;
    }
    
    return \%results;
});

$additional_io_timing{dna4_simd_input_scalar_output} = $dna4_simd_input_time;

print "\n";

# DNA2 Test 4: Sequential search without GIN index
print "=" x 70 . "\n";
print "DNA2 Test 4: Sequential search WITHOUT GIN index\n";
print "=" x 70 . "\n\n";

my %dna2_seq_search_results;
foreach my $simd (@simd_capabilities) {
    print "Testing DNA2 sequential search with force_simd_capability = $simd->{value} ($simd->{name})\n";
    
    # Set SIMD capability
    $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
    $dbh->do("SET kmersearch.kmer_size = $kmer_size");
    
    my ($seq_search_time, $seq_search_result) = measure_time("dna2_seq_search_$simd->{value}", sub {
        my $table_name = "simd_test_dna2_seq_" . $simd->{value};
        $table_name =~ s/\./_/g;
        
        create_base_table($table_name);
        
        # No GIN index creation - will force sequential scan
        
        # Store results for comparison
        my @all_results;
        my $batch_size = 10;  # Smaller batch size for sequential scan
        my $test_queries = 100;  # Fewer queries for sequential scan test
        
        # Test with subset of query sequences
        for (my $batch_start = 0; $batch_start < $test_queries && $batch_start < @query_sequences; $batch_start += $batch_size) {
            my $batch_end = $batch_start + $batch_size - 1;
            $batch_end = $test_queries - 1 if $batch_end >= $test_queries;
            $batch_end = $#query_sequences if $batch_end > $#query_sequences;
            
            for (my $i = $batch_start; $i <= $batch_end; $i++) {
                my $query = $query_sequences[$i];
                
                # Get matches with scores (sequential scan)
                my $sth = $dbh->prepare("
                    SELECT id, 
                           kmersearch_rawscore(seq, ?) as rawscore,
                           kmersearch_correctedscore(seq, ?) as correctedscore
                    FROM $table_name
                    WHERE seq =% ?
                    ORDER BY id
                ");
                $sth->execute($query, $query, $query);
                
                my @matches;
                while (my $row = $sth->fetchrow_hashref()) {
                    push @matches, {
                        id => $row->{id},
                        rawscore => $row->{rawscore},
                        correctedscore => $row->{correctedscore}
                    };
                }
                $sth->finish();
                
                push @all_results, {
                    query_idx => $i,
                    match_count => scalar(@matches),
                    matches => \@matches
                };
            }
            
            printf "    Processed %d/%d queries...\r", 
                   ($batch_end + 1), $test_queries;
        }
        print "\n";
        
        # Calculate summary statistics
        my $total_matches = 0;
        my $total_rawscore = 0;
        my $total_correctedscore = 0;
        my $queries_with_matches = 0;
        
        foreach my $result (@all_results) {
            $total_matches += $result->{match_count};
            $queries_with_matches++ if $result->{match_count} > 0;
            
            foreach my $match (@{$result->{matches}}) {
                $total_rawscore += $match->{rawscore};
                $total_correctedscore += $match->{correctedscore};
            }
        }
        
        # Clean up
        $dbh->do("DROP TABLE $table_name");
        
        return {
            query_count => scalar(@all_results),
            queries_with_matches => $queries_with_matches,
            total_matches => $total_matches,
            avg_matches_per_query => sprintf("%.2f", $total_matches / scalar(@all_results)),
            total_rawscore => $total_rawscore,
            total_correctedscore => $total_correctedscore,
            results_sample => [@all_results[0..9]]  # Keep first 10 for comparison
        };
    });
    
    $dna2_seq_search_results{$simd->{value}} = $seq_search_result;
    $timing_results{dna2_seq_search}{$simd->{value}} = $seq_search_time;
    printf "  Time: %.3fs, Total matches: %d, Queries with matches: %d/%d\n", 
           $seq_search_time, $seq_search_result->{total_matches}, 
           $seq_search_result->{queries_with_matches}, $seq_search_result->{query_count};
    printf "  Avg matches/query: %s, Total rawscore: %d, Total correctedscore: %d\n\n",
           $seq_search_result->{avg_matches_per_query}, 
           $seq_search_result->{total_rawscore}, $seq_search_result->{total_correctedscore};
}

# DNA4 Test 4: Sequential search without GIN index
print "=" x 70 . "\n";
print "DNA4 Test 4: Sequential search WITHOUT GIN index\n";
print "=" x 70 . "\n\n";

my %dna4_seq_search_results;
foreach my $simd (@simd_capabilities) {
    print "Testing DNA4 sequential search with force_simd_capability = $simd->{value} ($simd->{name})\n";
    
    # Set SIMD capability
    $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
    $dbh->do("SET kmersearch.kmer_size = $kmer_size");
    
    my ($seq_search_time, $seq_search_result) = measure_time("dna4_seq_search_$simd->{value}", sub {
        my $table_name = "simd_test_dna4_seq_" . $simd->{value};
        $table_name =~ s/\./_/g;
        
        create_base_table($table_name, 'DNA4');
        
        # No GIN index creation - will force sequential scan
        
        # Store results for comparison
        my @all_results;
        my $batch_size = 10;  # Smaller batch size for sequential scan
        my $test_queries = 100;  # Fewer queries for sequential scan test
        
        # Test with subset of query sequences
        for (my $batch_start = 0; $batch_start < $test_queries && $batch_start < @dna4_query_sequences; $batch_start += $batch_size) {
            my $batch_end = $batch_start + $batch_size - 1;
            $batch_end = $test_queries - 1 if $batch_end >= $test_queries;
            $batch_end = $#dna4_query_sequences if $batch_end > $#dna4_query_sequences;
            
            for (my $i = $batch_start; $i <= $batch_end; $i++) {
                my $query = $dna4_query_sequences[$i];
                
                # Get matches with scores (sequential scan)
                my $sth = $dbh->prepare("
                    SELECT id, 
                           kmersearch_rawscore(seq, ?) as rawscore,
                           kmersearch_correctedscore(seq, ?) as correctedscore
                    FROM $table_name
                    WHERE seq =% ?
                    ORDER BY id
                ");
                $sth->execute($query, $query, $query);
                
                my @matches;
                while (my $row = $sth->fetchrow_hashref()) {
                    push @matches, {
                        id => $row->{id},
                        rawscore => $row->{rawscore},
                        correctedscore => $row->{correctedscore}
                    };
                }
                $sth->finish();
                
                push @all_results, {
                    query_idx => $i,
                    match_count => scalar(@matches),
                    matches => \@matches
                };
            }
            
            printf "    Processed %d/%d queries...\r", 
                   ($batch_end + 1), $test_queries;
        }
        print "\n";
        
        # Calculate summary statistics
        my $total_matches = 0;
        my $total_rawscore = 0;
        my $total_correctedscore = 0;
        my $queries_with_matches = 0;
        
        foreach my $result (@all_results) {
            $total_matches += $result->{match_count};
            $queries_with_matches++ if $result->{match_count} > 0;
            
            foreach my $match (@{$result->{matches}}) {
                $total_rawscore += $match->{rawscore};
                $total_correctedscore += $match->{correctedscore};
            }
        }
        
        # Clean up
        $dbh->do("DROP TABLE $table_name");
        
        return {
            query_count => scalar(@all_results),
            queries_with_matches => $queries_with_matches,
            total_matches => $total_matches,
            avg_matches_per_query => sprintf("%.2f", $total_matches / scalar(@all_results)),
            total_rawscore => $total_rawscore,
            total_correctedscore => $total_correctedscore,
            results_sample => [@all_results[0..9]]  # Keep first 10 for comparison
        };
    });
    
    $dna4_seq_search_results{$simd->{value}} = $seq_search_result;
    $timing_results{dna4_seq_search}{$simd->{value}} = $seq_search_time;
    printf "  Time: %.3fs, Total matches: %d, Queries with matches: %d/%d\n", 
           $seq_search_time, $seq_search_result->{total_matches}, 
           $seq_search_result->{queries_with_matches}, $seq_search_result->{query_count};
    printf "  Avg matches/query: %s, Total rawscore: %d, Total correctedscore: %d\n\n",
           $seq_search_result->{avg_matches_per_query}, 
           $seq_search_result->{total_rawscore}, $seq_search_result->{total_correctedscore};
    
    # Test 5: Search with =% operator WITH GIN index
    print "\nTest 5: Search with =% operator WITH GIN index\n";
    my ($search_time, $search_result) = measure_time("search_scoring", sub {
        my $table_name = $gin_result->{table_name};
        
        # Store results for comparison
        my @all_results;
        my $batch_size = 100;  # Process queries in batches
        
        # Test with all 1000 query sequences
        for (my $batch_start = 0; $batch_start < @query_sequences; $batch_start += $batch_size) {
            my $batch_end = $batch_start + $batch_size - 1;
            $batch_end = $#query_sequences if $batch_end > $#query_sequences;
            
            for (my $i = $batch_start; $i <= $batch_end; $i++) {
                my $query = $query_sequences[$i];
                
                # Get matches with scores
                my $sth = $dbh->prepare("
                    SELECT id, 
                           kmersearch_rawscore(seq, ?) as rawscore,
                           kmersearch_correctedscore(seq, ?) as correctedscore
                    FROM $table_name
                    WHERE seq =% ?
                    ORDER BY id
                ");
                $sth->execute($query, $query, $query);
                
                my @matches;
                while (my $row = $sth->fetchrow_hashref()) {
                    push @matches, {
                        id => $row->{id},
                        rawscore => $row->{rawscore},
                        correctedscore => $row->{correctedscore}
                    };
                }
                $sth->finish();
                
                push @all_results, {
                    query_idx => $i,
                    match_count => scalar(@matches),
                    matches => \@matches
                };
            }
            
            printf "    Processed %d/%d queries...\r", 
                   ($batch_end + 1), scalar(@query_sequences);
        }
        print "\n";
        
        # Calculate summary statistics
        my $total_matches = 0;
        my $total_rawscore = 0;
        my $total_correctedscore = 0;
        my $queries_with_matches = 0;
        
        foreach my $result (@all_results) {
            $total_matches += $result->{match_count};
            $queries_with_matches++ if $result->{match_count} > 0;
            
            foreach my $match (@{$result->{matches}}) {
                $total_rawscore += $match->{rawscore};
                $total_correctedscore += $match->{correctedscore};
            }
        }
        
        # Clean up
        $dbh->do("DROP TABLE $table_name CASCADE");
        
        return {
            query_count => scalar(@all_results),
            queries_with_matches => $queries_with_matches,
            total_matches => $total_matches,
            avg_matches_per_query => sprintf("%.2f", $total_matches / scalar(@all_results)),
            total_rawscore => $total_rawscore,
            total_correctedscore => $total_correctedscore,
            results_sample => [@all_results[0..9]]  # Keep first 10 for comparison
        };
    });
    
    $timing_results{$simd->{value}}{search_scoring} = $search_time;
    $test_results{$simd->{value}}{search_scoring} = $search_result;
    printf "  Time: %.3fs, Total matches: %d, Queries with matches: %d/%d\n", 
           $search_time, $search_result->{total_matches}, 
           $search_result->{queries_with_matches}, $search_result->{query_count};
    printf "  Avg matches/query: %s, Total rawscore: %d, Total correctedscore: %d\n\n",
           $search_result->{avg_matches_per_query}, 
           $search_result->{total_rawscore}, $search_result->{total_correctedscore};
}

# DNA2 Test 5: K-mer frequency analysis
print "=" x 70 . "\n";
print "DNA2 Test 5: K-mer frequency analysis\n";
print "=" x 70 . "\n\n";

my %dna2_freq_results;
foreach my $simd (@simd_capabilities) {
    print "Testing DNA2 frequency analysis with force_simd_capability = $simd->{value} ($simd->{name})\n";
    
    # Set SIMD capability
    $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
    $dbh->do("SET kmersearch.kmer_size = $kmer_size");
    $dbh->do("SET kmersearch.max_appearance_rate = $max_appearance_rate");
    
    my ($freq_time, $freq_result) = measure_time("dna2_freq_$simd->{value}", sub {
        my $table_name = "simd_test_dna2_freq_" . $simd->{value};
        $table_name =~ s/\./_/g;
        
        create_base_table($table_name);
        
        # Perform analysis
        my $sth = $dbh->prepare("SELECT * FROM kmersearch_perform_highfreq_analysis(?, ?)");
        $sth->execute($table_name, 'seq');
        my $analysis = $sth->fetchrow_hashref();
        $sth->finish();
        
        # Get high-frequency k-mers for comparison
        $sth = $dbh->prepare("
            SELECT kmer2_as_uint::text, detection_reason 
            FROM kmersearch_highfreq_kmer 
            WHERE table_oid = ?::regclass::oid 
            ORDER BY kmer2_as_uint
            LIMIT 20
        ");
        $sth->execute($table_name);
        
        my @highfreq_kmers;
        while (my ($kmer, $reason) = $sth->fetchrow_array()) {
            push @highfreq_kmers, { kmer => $kmer, reason => $reason };
        }
        $sth->finish();
        
        # Get total count
        my ($total_highfreq) = $dbh->selectrow_array("
            SELECT COUNT(*) 
            FROM kmersearch_highfreq_kmer 
            WHERE table_oid = ?::regclass::oid
        ", undef, $table_name);
        
        # Clean up
        $dbh->do("SELECT kmersearch_undo_highfreq_analysis('$table_name', 'seq')");
        $dbh->do("DROP TABLE $table_name");
        
        return {
            total_rows => $analysis->{total_rows},
            highfreq_kmers_count => $analysis->{highfreq_kmers_count},
            highfreq_kmers_sample => \@highfreq_kmers,
            total_highfreq => $total_highfreq,
            parallel_workers => $analysis->{parallel_workers_used}
        };
    });
    
    $dna2_freq_results{$simd->{value}} = $freq_result;
    $timing_results{dna2_freq}{$simd->{value}} = $freq_time;
    printf "  Time: %.3fs, High-freq k-mers: %d, Workers: %d\n\n", 
           $freq_time, $freq_result->{highfreq_kmers_count}, $freq_result->{parallel_workers};
}

# DNA4 Test 5: K-mer frequency analysis
print "=" x 70 . "\n";
print "DNA4 Test 5: K-mer frequency analysis\n";
print "=" x 70 . "\n\n";

my %dna4_freq_results;
foreach my $simd (@simd_capabilities) {
    print "Testing DNA4 frequency analysis with force_simd_capability = $simd->{value} ($simd->{name})\n";
    
    # Set SIMD capability
    $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
    $dbh->do("SET kmersearch.kmer_size = $kmer_size");
    $dbh->do("SET kmersearch.max_appearance_rate = $max_appearance_rate");
    
    my ($freq_time, $freq_result) = measure_time("dna4_freq_$simd->{value}", sub {
        my $table_name = "simd_test_dna4_freq_" . $simd->{value};
        $table_name =~ s/\./_/g;
        
        create_base_table($table_name, 'DNA4');
        
        # Perform analysis
        my $sth = $dbh->prepare("SELECT * FROM kmersearch_perform_highfreq_analysis(?, ?)");
        $sth->execute($table_name, 'seq');
        my $analysis = $sth->fetchrow_hashref();
        $sth->finish();
        
        # Get high-frequency k-mers for comparison
        $sth = $dbh->prepare("
            SELECT kmer2_as_uint::text, detection_reason 
            FROM kmersearch_highfreq_kmer 
            WHERE table_oid = ?::regclass::oid 
            ORDER BY kmer2_as_uint
            LIMIT 20
        ");
        $sth->execute($table_name);
        
        my @highfreq_kmers;
        while (my ($kmer, $reason) = $sth->fetchrow_array()) {
            push @highfreq_kmers, { kmer => $kmer, reason => $reason };
        }
        $sth->finish();
        
        # Get total count
        my ($total_highfreq) = $dbh->selectrow_array("
            SELECT COUNT(*) 
            FROM kmersearch_highfreq_kmer 
            WHERE table_oid = ?::regclass::oid
        ", undef, $table_name);
        
        # Clean up
        $dbh->do("SELECT kmersearch_undo_highfreq_analysis('$table_name', 'seq')");
        $dbh->do("DROP TABLE $table_name");
        
        return {
            total_rows => $analysis->{total_rows},
            highfreq_kmers_count => $analysis->{highfreq_kmers_count},
            highfreq_kmers_sample => \@highfreq_kmers,
            total_highfreq => $total_highfreq,
            parallel_workers => $analysis->{parallel_workers_used}
        };
    });
    
    $dna4_freq_results{$simd->{value}} = $freq_result;
    $timing_results{dna4_freq}{$simd->{value}} = $freq_time;
    printf "  Time: %.3fs, High-freq k-mers: %d, Workers: %d\n\n", 
           $freq_time, $freq_result->{highfreq_kmers_count}, $freq_result->{parallel_workers};
}

# DNA2 Test 6: GIN index construction
print "=" x 70 . "\n";
print "DNA2 Test 6: GIN index construction\n";
print "=" x 70 . "\n\n";

my %dna2_gin_results;
my %dna2_gin_tables;  # Store table names for later use
foreach my $simd (@simd_capabilities) {
    print "Testing DNA2 GIN index with force_simd_capability = $simd->{value} ($simd->{name})\n";
    
    # Set SIMD capability
    $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
    $dbh->do("SET kmersearch.kmer_size = $kmer_size");
    $dbh->do("SET kmersearch.max_appearance_rate = $max_appearance_rate");
    $dbh->do("SET kmersearch.preclude_highfreq_kmer = true");
    
    my ($gin_time, $gin_result) = measure_time("dna2_gin_$simd->{value}", sub {
        my $table_name = "simd_test_dna2_gin_" . $simd->{value};
        $table_name =~ s/\./_/g;
        
        create_base_table($table_name);
        
        # Perform frequency analysis
        my ($analysis_rows) = $dbh->selectrow_array(
            "SELECT total_rows FROM kmersearch_perform_highfreq_analysis(?, ?)", 
            undef, $table_name, 'seq'
        );
        
        # Create GIN index
        my $index_name = $table_name . "_idx";
        my $gin_start = time();
        eval {
            local $SIG{ALRM} = sub { die "GIN index creation timeout\n" };
            alarm(60);  # 60 second timeout for GIN index creation
            $dbh->do("CREATE INDEX $index_name ON $table_name USING gin(seq kmersearch_dna2_gin_ops)");
            alarm(0);
        };
        my $gin_create_time = time() - $gin_start;
        if ($@) {
            print "  ERROR creating GIN index: $@\n";
            $dbh->do("DROP TABLE IF EXISTS $table_name CASCADE");
            return { error => $@ };
        }
        
        # Get index metadata
        my $sth = $dbh->prepare("
            SELECT highfreq_kmers_excluded 
            FROM kmersearch_gin_index_meta 
            WHERE index_oid = ?::regclass::oid
        ");
        $sth->execute($index_name);
        my ($excluded) = $sth->fetchrow_array();
        $sth->finish();
        
        return {
            table_name => $table_name,
            index_name => $index_name,
            highfreq_excluded => $excluded || 0,
            gin_create_time => $gin_create_time
        };
    });
    
    $dna2_gin_results{$simd->{value}} = $gin_result;
    $dna2_gin_tables{$simd->{value}} = $gin_result->{table_name};
    $timing_results{dna2_gin}{$simd->{value}} = $gin_time;
    printf "  Time: %.3fs (index creation: %.3fs), High-freq k-mers excluded: %d\n\n", 
           $gin_time, $gin_result->{gin_create_time}, $gin_result->{highfreq_excluded};
}

# DNA4 Test 6: GIN index construction
print "=" x 70 . "\n";
print "DNA4 Test 6: GIN index construction\n";
print "=" x 70 . "\n\n";

my %dna4_gin_results;
my %dna4_gin_tables;  # Store table names for later use
foreach my $simd (@simd_capabilities) {
    print "Testing DNA4 GIN index with force_simd_capability = $simd->{value} ($simd->{name})\n";
    
    # Set SIMD capability
    $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
    $dbh->do("SET kmersearch.kmer_size = $kmer_size");
    $dbh->do("SET kmersearch.max_appearance_rate = $max_appearance_rate");
    $dbh->do("SET kmersearch.preclude_highfreq_kmer = true");
    
    my ($gin_time, $gin_result) = measure_time("dna4_gin_$simd->{value}", sub {
        my $table_name = "simd_test_dna4_gin_" . $simd->{value};
        $table_name =~ s/\./_/g;
        
        create_base_table($table_name, 'DNA4');
        
        # Perform frequency analysis
        my ($analysis_rows) = $dbh->selectrow_array(
            "SELECT total_rows FROM kmersearch_perform_highfreq_analysis(?, ?)", 
            undef, $table_name, 'seq'
        );
        
        # Create GIN index
        my $index_name = $table_name . "_idx";
        my $gin_start = time();
        eval {
            local $SIG{ALRM} = sub { die "GIN index creation timeout\n" };
            alarm(60);  # 60 second timeout for GIN index creation
            $dbh->do("CREATE INDEX $index_name ON $table_name USING gin(seq kmersearch_dna4_gin_ops)");
            alarm(0);
        };
        my $gin_create_time = time() - $gin_start;
        if ($@) {
            print "  ERROR creating GIN index: $@\n";
            $dbh->do("DROP TABLE IF EXISTS $table_name CASCADE");
            return { error => $@ };
        }
        
        # Get index metadata
        my $sth = $dbh->prepare("
            SELECT highfreq_kmers_excluded 
            FROM kmersearch_gin_index_meta 
            WHERE index_oid = ?::regclass::oid
        ");
        $sth->execute($index_name);
        my ($excluded) = $sth->fetchrow_array();
        $sth->finish();
        
        return {
            table_name => $table_name,
            index_name => $index_name,
            highfreq_excluded => $excluded || 0,
            gin_create_time => $gin_create_time
        };
    });
    
    $dna4_gin_results{$simd->{value}} = $gin_result;
    $dna4_gin_tables{$simd->{value}} = $gin_result->{table_name};
    $timing_results{dna4_gin}{$simd->{value}} = $gin_time;
    printf "  Time: %.3fs (index creation: %.3fs), High-freq k-mers excluded: %d\n\n", 
           $gin_time, $gin_result->{gin_create_time}, $gin_result->{highfreq_excluded};
}

# DNA2 Test 7: GIN index search
print "=" x 70 . "\n";
print "DNA2 Test 7: GIN index search\n";
print "=" x 70 . "\n\n";

my %dna2_gin_search_results;
foreach my $simd (@simd_capabilities) {
    print "Testing DNA2 GIN search with force_simd_capability = $simd->{value} ($simd->{name})\n";
    
    # Set SIMD capability
    $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
    $dbh->do("SET kmersearch.kmer_size = $kmer_size");
    
    my ($search_time, $search_result) = measure_time("dna2_gin_search_$simd->{value}", sub {
        my $table_name = $dna2_gin_tables{$simd->{value}};
        
        # Store results for comparison
        my @all_results;
        my $batch_size = 100;  # Process queries in batches
        
        # Test with all 1000 query sequences
        for (my $batch_start = 0; $batch_start < @query_sequences; $batch_start += $batch_size) {
            my $batch_end = $batch_start + $batch_size - 1;
            $batch_end = $#query_sequences if $batch_end > $#query_sequences;
            
            for (my $i = $batch_start; $i <= $batch_end; $i++) {
                my $query = $query_sequences[$i];
                
                # Get matches with scores
                my $sth = $dbh->prepare("
                    SELECT id, 
                           kmersearch_rawscore(seq, ?) as rawscore,
                           kmersearch_correctedscore(seq, ?) as correctedscore
                    FROM $table_name
                    WHERE seq =% ?
                    ORDER BY id
                ");
                $sth->execute($query, $query, $query);
                
                my @matches;
                while (my $row = $sth->fetchrow_hashref()) {
                    push @matches, {
                        id => $row->{id},
                        rawscore => $row->{rawscore},
                        correctedscore => $row->{correctedscore}
                    };
                }
                $sth->finish();
                
                push @all_results, {
                    query_idx => $i,
                    match_count => scalar(@matches),
                    matches => \@matches
                };
            }
            
            printf "    Processed %d/%d queries...\r", 
                   ($batch_end + 1), scalar(@query_sequences);
        }
        print "\n";
        
        # Calculate summary statistics
        my $total_matches = 0;
        my $total_rawscore = 0;
        my $total_correctedscore = 0;
        my $queries_with_matches = 0;
        
        foreach my $result (@all_results) {
            $total_matches += $result->{match_count};
            $queries_with_matches++ if $result->{match_count} > 0;
            
            foreach my $match (@{$result->{matches}}) {
                $total_rawscore += $match->{rawscore};
                $total_correctedscore += $match->{correctedscore};
            }
        }
        
        # Clean up
        $dbh->do("DROP TABLE $table_name CASCADE");
        
        return {
            query_count => scalar(@all_results),
            queries_with_matches => $queries_with_matches,
            total_matches => $total_matches,
            avg_matches_per_query => sprintf("%.2f", $total_matches / scalar(@all_results)),
            total_rawscore => $total_rawscore,
            total_correctedscore => $total_correctedscore,
            results_sample => [@all_results[0..9]]  # Keep first 10 for comparison
        };
    });
    
    $dna2_gin_search_results{$simd->{value}} = $search_result;
    $timing_results{dna2_gin_search}{$simd->{value}} = $search_time;
    printf "  Time: %.3fs, Total matches: %d, Queries with matches: %d/%d\n", 
           $search_time, $search_result->{total_matches}, 
           $search_result->{queries_with_matches}, $search_result->{query_count};
    printf "  Avg matches/query: %s, Total rawscore: %d, Total correctedscore: %d\n\n",
           $search_result->{avg_matches_per_query}, 
           $search_result->{total_rawscore}, $search_result->{total_correctedscore};
}

# DNA4 Test 7: GIN index search
print "=" x 70 . "\n";
print "DNA4 Test 7: GIN index search\n";
print "=" x 70 . "\n\n";

my %dna4_gin_search_results;
foreach my $simd (@simd_capabilities) {
    print "Testing DNA4 GIN search with force_simd_capability = $simd->{value} ($simd->{name})\n";
    
    # Set SIMD capability
    $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
    $dbh->do("SET kmersearch.kmer_size = $kmer_size");
    
    my ($search_time, $search_result) = measure_time("dna4_gin_search_$simd->{value}", sub {
        my $table_name = $dna4_gin_tables{$simd->{value}};
        
        # Store results for comparison
        my @all_results;
        my $batch_size = 100;  # Process queries in batches
        
        # Test with all 1000 query sequences
        for (my $batch_start = 0; $batch_start < @dna4_query_sequences; $batch_start += $batch_size) {
            my $batch_end = $batch_start + $batch_size - 1;
            $batch_end = $#dna4_query_sequences if $batch_end > $#dna4_query_sequences;
            
            for (my $i = $batch_start; $i <= $batch_end; $i++) {
                my $query = $dna4_query_sequences[$i];
                
                # Get matches with scores
                my $sth = $dbh->prepare("
                    SELECT id, 
                           kmersearch_rawscore(seq, ?) as rawscore,
                           kmersearch_correctedscore(seq, ?) as correctedscore
                    FROM $table_name
                    WHERE seq =% ?
                    ORDER BY id
                ");
                $sth->execute($query, $query, $query);
                
                my @matches;
                while (my $row = $sth->fetchrow_hashref()) {
                    push @matches, {
                        id => $row->{id},
                        rawscore => $row->{rawscore},
                        correctedscore => $row->{correctedscore}
                    };
                }
                $sth->finish();
                
                push @all_results, {
                    query_idx => $i,
                    match_count => scalar(@matches),
                    matches => \@matches
                };
            }
            
            printf "    Processed %d/%d queries...\r", 
                   ($batch_end + 1), scalar(@dna4_query_sequences);
        }
        print "\n";
        
        # Calculate summary statistics
        my $total_matches = 0;
        my $total_rawscore = 0;
        my $total_correctedscore = 0;
        my $queries_with_matches = 0;
        
        foreach my $result (@all_results) {
            $total_matches += $result->{match_count};
            $queries_with_matches++ if $result->{match_count} > 0;
            
            foreach my $match (@{$result->{matches}}) {
                $total_rawscore += $match->{rawscore};
                $total_correctedscore += $match->{correctedscore};
            }
        }
        
        # Clean up
        $dbh->do("DROP TABLE $table_name CASCADE");
        
        return {
            query_count => scalar(@all_results),
            queries_with_matches => $queries_with_matches,
            total_matches => $total_matches,
            avg_matches_per_query => sprintf("%.2f", $total_matches / scalar(@all_results)),
            total_rawscore => $total_rawscore,
            total_correctedscore => $total_correctedscore,
            results_sample => [@all_results[0..9]]  # Keep first 10 for comparison
        };
    });
    
    $dna4_gin_search_results{$simd->{value}} = $search_result;
    $timing_results{dna4_gin_search}{$simd->{value}} = $search_time;
    printf "  Time: %.3fs, Total matches: %d, Queries with matches: %d/%d\n", 
           $search_time, $search_result->{total_matches}, 
           $search_result->{queries_with_matches}, $search_result->{query_count};
    printf "  Avg matches/query: %s, Total rawscore: %d, Total correctedscore: %d\n\n",
           $search_result->{avg_matches_per_query}, 
           $search_result->{total_rawscore}, $search_result->{total_correctedscore};
}

# Compare results across SIMD capabilities
print "=" x 70 . "\n";
print "Result Comparison\n";
print "=" x 70 . "\n\n";

# Compare DNA2 sequential search results
print "DNA2 Sequential Search Results Comparison:\n";
my $base_dna2_seq = $dna2_seq_search_results{$simd_capabilities[0]->{value}};
my $dna2_seq_consistent = 1;

foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
    my $seq_search = $dna2_seq_search_results{$simd->{value}};
    
    if ($seq_search->{total_matches} != $base_dna2_seq->{total_matches} ||
        $seq_search->{total_rawscore} != $base_dna2_seq->{total_rawscore} ||
        $seq_search->{total_correctedscore} != $base_dna2_seq->{total_correctedscore}) {
        printf "  WARNING: Different DNA2 sequential search results for %s:\n", $simd->{name};
        printf "    Matches: %d vs %d\n", $seq_search->{total_matches}, $base_dna2_seq->{total_matches};
        printf "    Rawscore: %d vs %d\n", $seq_search->{total_rawscore}, $base_dna2_seq->{total_rawscore};
        printf "    Correctedscore: %d vs %d\n", $seq_search->{total_correctedscore}, $base_dna2_seq->{total_correctedscore};
        $dna2_seq_consistent = 0;
    }
}

if ($dna2_seq_consistent) {
    print "  ✓ DNA2 Sequential search results are consistent across all SIMD capabilities\n";
}

# Compare DNA4 sequential search results
print "\nDNA4 Sequential Search Results Comparison:\n";
my $base_dna4_seq = $dna4_seq_search_results{$simd_capabilities[0]->{value}};
my $dna4_seq_consistent = 1;

foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
    my $seq_search = $dna4_seq_search_results{$simd->{value}};
    
    if ($seq_search->{total_matches} != $base_dna4_seq->{total_matches} ||
        $seq_search->{total_rawscore} != $base_dna4_seq->{total_rawscore} ||
        $seq_search->{total_correctedscore} != $base_dna4_seq->{total_correctedscore}) {
        printf "  WARNING: Different DNA4 sequential search results for %s:\n", $simd->{name};
        printf "    Matches: %d vs %d\n", $seq_search->{total_matches}, $base_dna4_seq->{total_matches};
        printf "    Rawscore: %d vs %d\n", $seq_search->{total_rawscore}, $base_dna4_seq->{total_rawscore};
        printf "    Correctedscore: %d vs %d\n", $seq_search->{total_correctedscore}, $base_dna4_seq->{total_correctedscore};
        $dna4_seq_consistent = 0;
    }
}

if ($dna4_seq_consistent) {
    print "  ✓ DNA4 Sequential search results are consistent across all SIMD capabilities\n";
}

# Compare DNA2 frequency analysis results
print "\nDNA2 Frequency Analysis Comparison:\n";
my $base_dna2_freq = $dna2_freq_results{$simd_capabilities[0]->{value}};
my $dna2_freq_consistent = 1;

foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
    my $freq = $dna2_freq_results{$simd->{value}};
    
    if ($freq->{highfreq_kmers_count} != $base_dna2_freq->{highfreq_kmers_count}) {
        printf "  WARNING: Different DNA2 high-freq k-mer count for %s: %d vs %d\n",
               $simd->{name}, $freq->{highfreq_kmers_count}, $base_dna2_freq->{highfreq_kmers_count};
        $dna2_freq_consistent = 0;
    }
}

if ($dna2_freq_consistent) {
    print "  ✓ DNA2 High-frequency k-mer analysis results are consistent across all SIMD capabilities\n";
}

# Compare DNA4 frequency analysis results
print "\nDNA4 Frequency Analysis Comparison:\n";
my $base_dna4_freq = $dna4_freq_results{$simd_capabilities[0]->{value}};
my $dna4_freq_consistent = 1;

foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
    my $freq = $dna4_freq_results{$simd->{value}};
    
    if ($freq->{highfreq_kmers_count} != $base_dna4_freq->{highfreq_kmers_count}) {
        printf "  WARNING: Different DNA4 high-freq k-mer count for %s: %d vs %d\n",
               $simd->{name}, $freq->{highfreq_kmers_count}, $base_dna4_freq->{highfreq_kmers_count};
        $dna4_freq_consistent = 0;
    }
}

if ($dna4_freq_consistent) {
    print "  ✓ DNA4 High-frequency k-mer analysis results are consistent across all SIMD capabilities\n";
}

# Compare DNA2 GIN index search results
print "\nDNA2 GIN Index Search Results Comparison:\n";
my $base_dna2_gin_search = $dna2_gin_search_results{$simd_capabilities[0]->{value}};
my $dna2_gin_search_consistent = 1;

foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
    my $search = $dna2_gin_search_results{$simd->{value}};
    
    if ($search->{total_matches} != $base_dna2_gin_search->{total_matches} ||
        $search->{total_rawscore} != $base_dna2_gin_search->{total_rawscore} ||
        $search->{total_correctedscore} != $base_dna2_gin_search->{total_correctedscore}) {
        printf "  WARNING: Different DNA2 GIN search results for %s:\n", $simd->{name};
        printf "    Matches: %d vs %d\n", $search->{total_matches}, $base_dna2_gin_search->{total_matches};
        printf "    Rawscore: %d vs %d\n", $search->{total_rawscore}, $base_dna2_gin_search->{total_rawscore};
        printf "    Correctedscore: %d vs %d\n", $search->{total_correctedscore}, $base_dna2_gin_search->{total_correctedscore};
        $dna2_gin_search_consistent = 0;
    }
}

if ($dna2_gin_search_consistent) {
    print "  ✓ DNA2 GIN index search results are consistent across all SIMD capabilities\n";
}

# Compare DNA4 GIN index search results
print "\nDNA4 GIN Index Search Results Comparison:\n";
my $base_dna4_gin_search = $dna4_gin_search_results{$simd_capabilities[0]->{value}};
my $dna4_gin_search_consistent = 1;

foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
    my $search = $dna4_gin_search_results{$simd->{value}};
    
    if ($search->{total_matches} != $base_dna4_gin_search->{total_matches} ||
        $search->{total_rawscore} != $base_dna4_gin_search->{total_rawscore} ||
        $search->{total_correctedscore} != $base_dna4_gin_search->{total_correctedscore}) {
        printf "  WARNING: Different DNA4 GIN search results for %s:\n", $simd->{name};
        printf "    Matches: %d vs %d\n", $search->{total_matches}, $base_dna4_gin_search->{total_matches};
        printf "    Rawscore: %d vs %d\n", $search->{total_rawscore}, $base_dna4_gin_search->{total_rawscore};
        printf "    Correctedscore: %d vs %d\n", $search->{total_correctedscore}, $base_dna4_gin_search->{total_correctedscore};
        $dna4_gin_search_consistent = 0;
    }
}

if ($dna4_gin_search_consistent) {
    print "  ✓ DNA4 GIN index search results are consistent across all SIMD capabilities\n";
}

# Performance Summary
print "\n" . "=" x 100 . "\n";
print "Performance Summary\n";
print "=" x 100 . "\n\n";

# Additional I/O Tests Performance
print "Additional I/O Tests:\n";
print "-" x 80 . "\n";
printf "%-30s %12s\n", "Test", "Time (s)";
print "-" x 80 . "\n";
printf "%-30s %12.3f\n", "DNA2 Scalar I/O", $additional_io_timing{dna2_scalar_io};
printf "%-30s %12.3f\n", "DNA4 Scalar I/O", $additional_io_timing{dna4_scalar_io};
printf "%-30s %12.3f\n", "DNA2 Scalar Input/SIMD Output", $additional_io_timing{dna2_scalar_input_simd_output};
printf "%-30s %12.3f\n", "DNA4 Scalar Input/SIMD Output", $additional_io_timing{dna4_scalar_input_simd_output};
printf "%-30s %12.3f\n", "DNA2 SIMD Input/Scalar Output", $additional_io_timing{dna2_simd_input_scalar_output};
printf "%-30s %12.3f\n", "DNA4 SIMD Input/Scalar Output", $additional_io_timing{dna4_simd_input_scalar_output};

# Main Test Performance Summary by SIMD capability
print "\nDNA2 Test Performance by SIMD Capability:\n";
print "-" x 90 . "\n";
printf "%-15s %12s %12s %12s %12s %12s\n", 
       "SIMD Capability", "Seq Search(s)", "Freq Anal(s)", "GIN Index(s)", "GIN Search(s)", "Total (s)";
print "-" x 90 . "\n";

foreach my $simd (@simd_capabilities) {
    my $total = $timing_results{dna2_seq_search}{$simd->{value}} +
                $timing_results{dna2_freq}{$simd->{value}} +
                $timing_results{dna2_gin}{$simd->{value}} +
                $timing_results{dna2_gin_search}{$simd->{value}};
    
    printf "%-15s %12.3f %12.3f %12.3f %12.3f %12.3f\n",
        $simd->{name},
        $timing_results{dna2_seq_search}{$simd->{value}},
        $timing_results{dna2_freq}{$simd->{value}},
        $timing_results{dna2_gin}{$simd->{value}},
        $timing_results{dna2_gin_search}{$simd->{value}},
        $total;
}

print "\nDNA4 Test Performance by SIMD Capability:\n";
print "-" x 90 . "\n";
printf "%-15s %12s %12s %12s %12s %12s\n", 
       "SIMD Capability", "Seq Search(s)", "Freq Anal(s)", "GIN Index(s)", "GIN Search(s)", "Total (s)";
print "-" x 90 . "\n";

foreach my $simd (@simd_capabilities) {
    my $total = $timing_results{dna4_seq_search}{$simd->{value}} +
                $timing_results{dna4_freq}{$simd->{value}} +
                $timing_results{dna4_gin}{$simd->{value}} +
                $timing_results{dna4_gin_search}{$simd->{value}};
    
    printf "%-15s %12.3f %12.3f %12.3f %12.3f %12.3f\n",
        $simd->{name},
        $timing_results{dna4_seq_search}{$simd->{value}},
        $timing_results{dna4_freq}{$simd->{value}},
        $timing_results{dna4_gin}{$simd->{value}},
        $timing_results{dna4_gin_search}{$simd->{value}},
        $total;
}

# Calculate speedup
if (@simd_capabilities > 1) {
    print "\nDNA2 Speedup relative to 'None':\n";
    print "-" x 90 . "\n";
    printf "%-15s %12s %12s %12s %12s %12s\n",
           "SIMD Capability", "Seq Search", "Freq Anal", "GIN Index", "GIN Search", "Total";
    print "-" x 90 . "\n";
    
    my $base_value = $simd_capabilities[0]->{value};
    my $base_dna2_total = $timing_results{dna2_seq_search}{$base_value} +
                          $timing_results{dna2_freq}{$base_value} +
                          $timing_results{dna2_gin}{$base_value} +
                          $timing_results{dna2_gin_search}{$base_value};
    
    foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
        my $curr_total = $timing_results{dna2_seq_search}{$simd->{value}} +
                         $timing_results{dna2_freq}{$simd->{value}} +
                         $timing_results{dna2_gin}{$simd->{value}} +
                         $timing_results{dna2_gin_search}{$simd->{value}};
        
        printf "%-15s %12.2fx %12.2fx %12.2fx %12.2fx %12.2fx\n",
            $simd->{name},
            $timing_results{dna2_seq_search}{$base_value} / $timing_results{dna2_seq_search}{$simd->{value}},
            $timing_results{dna2_freq}{$base_value} / $timing_results{dna2_freq}{$simd->{value}},
            $timing_results{dna2_gin}{$base_value} / $timing_results{dna2_gin}{$simd->{value}},
            $timing_results{dna2_gin_search}{$base_value} / $timing_results{dna2_gin_search}{$simd->{value}},
            $base_dna2_total / $curr_total;
    }
    
    print "\nDNA4 Speedup relative to 'None':\n";
    print "-" x 90 . "\n";
    printf "%-15s %12s %12s %12s %12s %12s\n",
           "SIMD Capability", "Seq Search", "Freq Anal", "GIN Index", "GIN Search", "Total";
    print "-" x 90 . "\n";
    
    my $base_dna4_total = $timing_results{dna4_seq_search}{$base_value} +
                          $timing_results{dna4_freq}{$base_value} +
                          $timing_results{dna4_gin}{$base_value} +
                          $timing_results{dna4_gin_search}{$base_value};
    
    foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
        my $curr_total = $timing_results{dna4_seq_search}{$simd->{value}} +
                         $timing_results{dna4_freq}{$simd->{value}} +
                         $timing_results{dna4_gin}{$simd->{value}} +
                         $timing_results{dna4_gin_search}{$simd->{value}};
        
        printf "%-15s %12.2fx %12.2fx %12.2fx %12.2fx %12.2fx\n",
            $simd->{name},
            $timing_results{dna4_seq_search}{$base_value} / $timing_results{dna4_seq_search}{$simd->{value}},
            $timing_results{dna4_freq}{$base_value} / $timing_results{dna4_freq}{$simd->{value}},
            $timing_results{dna4_gin}{$base_value} / $timing_results{dna4_gin}{$simd->{value}},
            $timing_results{dna4_gin_search}{$base_value} / $timing_results{dna4_gin_search}{$simd->{value}},
            $base_dna4_total / $curr_total;
    }
}

# Reset to auto-detect
print "\nResetting to auto-detected SIMD capability...\n";
$dbh->do("SET kmersearch.force_simd_capability = -1");
my ($final) = $dbh->selectrow_array("SELECT kmersearch_simd_capability()");
print "Final SIMD capability: $final\n";

print "\nAll tests completed successfully!\n";

# Disconnect
$dbh->disconnect();
print "=" x 70 . "\n";
print "Additional Input/Output Tests (DNA4)\n";
print "=" x 70 . "\n\n";

# DNA4 Test 1: Scalar input/Scalar output
print "DNA4 Test 1: Scalar input (force_simd_capability=0) / Scalar output (force_simd_capability=0)\n";
print "  Testing with " . scalar(@dna4_sequences) . " DNA4 sequences...\n";
my ($dna4_scalar_io_time, $dna4_scalar_io_result) = measure_time("dna4_scalar_io", sub {
    my $table_name = "simd_test_dna4_scalar_io";
    
    # Set to scalar mode for input
    $dbh->do("SET kmersearch.force_simd_capability = 0");
    
    # Create table and insert data
    create_base_table($table_name, 'DNA4');
    
    # Verify with scalar output
    $dbh->do("SET kmersearch.force_simd_capability = 0");
    
    my $sth = $dbh->prepare("SELECT id, seq::text FROM $table_name ORDER BY id");
    $sth->execute();
    
    my $match_count = 0;
    my $row_count = 0;
    while (my ($id, $seq) = $sth->fetchrow_array()) {
        if ($dna4_sequences[$id-1] eq $seq) {
            $match_count++;
        } else {
            print "  Mismatch at ID $id: expected '$dna4_sequences[$id-1]', got '$seq'\n" if $row_count < 5;
        }
        $row_count++;
    }
    $sth->finish();
    
    # Clean up
    $dbh->do("DROP TABLE $table_name");
    
    return {
        total_rows => $row_count,
        matched => $match_count,
        success => ($match_count == scalar(@dna4_sequences))
    };
});

$dna4_additional_io_timing{scalar_io} = $dna4_scalar_io_time;
printf "  Time: %.3fs, Result: %d/%d sequences match (%.1f%%)\n", 
       $dna4_scalar_io_time, $dna4_scalar_io_result->{matched}, $dna4_scalar_io_result->{total_rows},
       ($dna4_scalar_io_result->{matched} / $dna4_scalar_io_result->{total_rows}) * 100;

# DNA4 Test 2: Scalar input/SIMD output (all SIMD values)
print "\nDNA4 Test 2: Scalar input (force_simd_capability=0) / SIMD output (all force_simd_capability>0)\n";
print "  Testing " . scalar(@simd_capabilities) . " SIMD capabilities...\n";
my ($dna4_scalar_input_time, $dna4_scalar_input_result) = measure_time("dna4_scalar_input_simd_output", sub {
    my $table_name = "simd_test_dna4_scalar_input";
    my %results;
    
    # Set to scalar mode for input
    $dbh->do("SET kmersearch.force_simd_capability = 0");
    
    # Create table and insert data
    create_base_table($table_name, 'DNA4');
    
    # Test with each SIMD capability for output
    foreach my $simd (@simd_capabilities) {
        next if $simd->{value} == 0;  # Skip scalar
        
        my $test_start = time();
        $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
        
        my $sth = $dbh->prepare("SELECT id, seq::text FROM $table_name ORDER BY id");
        $sth->execute();
        
        my $match_count = 0;
        my $row_count = 0;
        while (my ($id, $seq) = $sth->fetchrow_array()) {
            if ($dna4_sequences[$id-1] eq $seq) {
                $match_count++;
            } else {
                print "  Mismatch with $simd->{name} output at ID $id\n" if $row_count < 5;
            }
            $row_count++;
        }
        $sth->finish();
        my $test_time = time() - $test_start;
        
        $results{$simd->{name}} = {
            time => $test_time,
            matched => $match_count,
            total => $row_count,
            success_rate => ($match_count / $row_count) * 100
        };
        
        printf "  $simd->{name} output: %d/%d sequences match (%.1f%%) - Time: %.3fs\n", 
               $match_count, $row_count, ($match_count / $row_count) * 100, $test_time;
    }
    
    # Clean up
    $dbh->do("DROP TABLE $table_name");
    
    return \%results;
});

$dna4_additional_io_timing{scalar_input_simd_output} = $dna4_scalar_input_time;

# DNA4 Test 3: SIMD input/Scalar output (all SIMD values)
print "\nDNA4 Test 3: SIMD input (all force_simd_capability>0) / Scalar output (force_simd_capability=0)\n";
print "  Testing " . scalar(@simd_capabilities) . " SIMD capabilities for input...\n";
my ($dna4_simd_input_time, $dna4_simd_input_result) = measure_time("dna4_simd_input_scalar_output", sub {
    my %results;
    
    foreach my $simd (@simd_capabilities) {
        next if $simd->{value} == 0;  # Skip scalar
        
        my $test_start = time();
        my $table_name = "simd_test_dna4_" . $simd->{name} . "_input";
        $table_name =~ s/[^a-zA-Z0-9_]/_/g;  # Replace non-alphanumeric chars
        
        # Set SIMD mode for input
        $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
        
        # Create table and insert data
        create_base_table($table_name, 'DNA4');
        
        # Set to scalar mode for output
        $dbh->do("SET kmersearch.force_simd_capability = 0");
        
        my $sth = $dbh->prepare("SELECT id, seq::text FROM $table_name ORDER BY id");
        $sth->execute();
        
        my $match_count = 0;
        my $row_count = 0;
        while (my ($id, $seq) = $sth->fetchrow_array()) {
            if ($dna4_sequences[$id-1] eq $seq) {
                $match_count++;
            } else {
                print "  Mismatch with $simd->{name} input at ID $id\n" if $row_count < 5;
            }
            $row_count++;
        }
        $sth->finish();
        
        # Clean up
        $dbh->do("DROP TABLE $table_name");
        
        my $test_time = time() - $test_start;
        
        $results{$simd->{name}} = {
            time => $test_time,
            matched => $match_count,
            total => $row_count,
            success_rate => ($match_count / $row_count) * 100
        };
        
        printf "  $simd->{name} input: %d/%d sequences match (%.1f%%) - Time: %.3fs\n", 
               $match_count, $row_count, ($match_count / $row_count) * 100, $test_time;
    }
    
    return \%results;
});

$dna4_additional_io_timing{simd_input_scalar_output} = $dna4_simd_input_time;

print "\n";

# Main DNA4 Tests
print "=" x 70 . "\n";
print "Main DNA4 Tests (per SIMD capability)\n";
print "=" x 70 . "\n\n";

# Test each SIMD capability with DNA4
foreach my $simd (@simd_capabilities) {
    print "=" x 70 . "\n";
    print "DNA4 Testing with kmersearch.force_simd_capability = '$simd->{value}'\n";
    print "=" x 70 . "\n";
    
    # Set SIMD capability
    $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
    $dbh->do("SET kmersearch.kmer_size = $kmer_size");
    
    # Verify setting
    my ($current) = $dbh->selectrow_array("SELECT kmersearch_simd_capability()");
    print "Current SIMD capability: $current\n\n";
    
    # DNA4 Test 1: Input/Output verification
    print "DNA4 Test 1: Input/Output verification\n";
    my ($io_time, $io_result) = measure_time("dna4_input_output", sub {
        my $table_name = "simd_test_dna4_io_" . $simd->{value};
        $table_name =~ s/\./_/g;  # Replace dots with underscores
        
        create_base_table($table_name, 'DNA4');
        
        # Verify data integrity
        my $sth = $dbh->prepare("SELECT id, seq::text FROM $table_name ORDER BY id");
        $sth->execute();
        
        my $match_count = 0;
        my $row_count = 0;
        while (my ($id, $seq) = $sth->fetchrow_array()) {
            if ($dna4_sequences[$id-1] eq $seq) {
                $match_count++;
            } else {
                print "  Mismatch at ID $id: expected '$dna4_sequences[$id-1]', got '$seq'\n" if $row_count < 5;
            }
            $row_count++;
        }
        $sth->finish();
        
        # Clean up
        $dbh->do("DROP TABLE $table_name");
        
        return {
            total_rows => $row_count,
            matched => $match_count,
            success => ($match_count == scalar(@dna4_sequences))
        };
    });
    
    $dna4_timing_results{$simd->{value}}{input_output} = $io_time;
    $dna4_test_results{$simd->{value}}{input_output} = $io_result;
    printf "  Time: %.3fs, Matched: %d/%d (%.1f%%)\n", 
           $io_time, $io_result->{matched}, $io_result->{total_rows},
           ($io_result->{matched} / $io_result->{total_rows}) * 100;
    
    # DNA4 Test 2: K-mer frequency analysis
    print "\nDNA4 Test 2: K-mer frequency analysis (max_appearance_rate=$max_appearance_rate)\n";
    my ($freq_time, $freq_result) = measure_time("dna4_frequency_analysis", sub {
        my $table_name = "simd_test_dna4_freq_" . $simd->{value};
        $table_name =~ s/\./_/g;
        
        create_base_table($table_name, 'DNA4');
        
        # Set parameters
        $dbh->do("SET kmersearch.max_appearance_rate = $max_appearance_rate");
        
        # Perform analysis
        my $sth = $dbh->prepare("SELECT * FROM kmersearch_perform_highfreq_analysis(?, ?)");
        $sth->execute($table_name, 'seq');
        my $analysis = $sth->fetchrow_hashref();
        $sth->finish();
        
        # Get high-frequency k-mers for comparison
        $sth = $dbh->prepare("
            SELECT kmer2_as_uint::text, detection_reason 
            FROM kmersearch_highfreq_kmer 
            WHERE table_oid = ?::regclass::oid 
            ORDER BY kmer2_as_uint
            LIMIT 20
        ");
        $sth->execute($table_name);
        
        my @highfreq_kmers;
        while (my ($kmer, $reason) = $sth->fetchrow_array()) {
            push @highfreq_kmers, { kmer => $kmer, reason => $reason };
        }
        $sth->finish();
        
        # Get total count
        my ($total_highfreq) = $dbh->selectrow_array("
            SELECT COUNT(*) 
            FROM kmersearch_highfreq_kmer 
            WHERE table_oid = ?::regclass::oid
        ", undef, $table_name);
        
        # Clean up
        $dbh->do("SELECT kmersearch_undo_highfreq_analysis('$table_name', 'seq')");
        $dbh->do("DROP TABLE $table_name");
        
        return {
            total_rows => $analysis->{total_rows},
            highfreq_kmers_count => $analysis->{highfreq_kmers_count},
            highfreq_kmers_sample => \@highfreq_kmers,
            total_highfreq => $total_highfreq,
            parallel_workers => $analysis->{parallel_workers_used}
        };
    });
    
    $dna4_timing_results{$simd->{value}}{frequency_analysis} = $freq_time;
    $dna4_test_results{$simd->{value}}{frequency_analysis} = $freq_result;
    printf "  Time: %.3fs, High-freq k-mers: %d, Workers: %d\n", 
           $freq_time, $freq_result->{highfreq_kmers_count}, $freq_result->{parallel_workers};
    
    # DNA4 Test 3: GIN index construction with high-frequency k-mer exclusion
    print "\nDNA4 Test 3: GIN index construction with high-frequency k-mer exclusion\n";
    my ($gin_time, $gin_result) = measure_time("dna4_gin_index", sub {
        my $table_name = "simd_test_dna4_gin_" . $simd->{value};
        $table_name =~ s/\./_/g;
        
        create_base_table($table_name, 'DNA4');
        
        # Set parameters
        $dbh->do("SET kmersearch.max_appearance_rate = $max_appearance_rate");
        $dbh->do("SET kmersearch.preclude_highfreq_kmer = true");
        
        # Perform frequency analysis
        my ($analysis_rows) = $dbh->selectrow_array(
            "SELECT total_rows FROM kmersearch_perform_highfreq_analysis(?, ?)", 
            undef, $table_name, 'seq'
        );
        
        # Create GIN index
        my $index_name = $table_name . "_idx";
        print "  Creating GIN index $index_name...\n";
        my $gin_start = time();
        eval {
            local $SIG{ALRM} = sub { die "GIN index creation timeout\n" };
            alarm(60);  # 60 second timeout for GIN index creation
            $dbh->do("CREATE INDEX $index_name ON $table_name USING gin(seq kmersearch_dna4_gin_ops)");
            alarm(0);
        };
        my $gin_create_time = time() - $gin_start;
        if ($@) {
            print "  ERROR creating GIN index: $@\n";
            $dbh->do("DROP TABLE IF EXISTS $table_name CASCADE");
            return { error => $@ };
        }
        print "  GIN index created in ${gin_create_time}s\n";
        
        # Get index metadata
        my $sth = $dbh->prepare("
            SELECT highfreq_kmers_excluded 
            FROM kmersearch_gin_index_meta 
            WHERE index_oid = ?::regclass::oid
        ");
        $sth->execute($index_name);
        my ($excluded) = $sth->fetchrow_array();
        $sth->finish();
        
        return {
            table_name => $table_name,
            index_name => $index_name,
            highfreq_excluded => $excluded || 0,
            gin_create_time => $gin_create_time
        };
    });
    
    $dna4_timing_results{$simd->{value}}{gin_index} = $gin_time;
    $dna4_test_results{$simd->{value}}{gin_index} = $gin_result;
    printf "  Time: %.3fs (index creation: %.3fs), High-freq k-mers excluded: %d\n", 
           $gin_time, $gin_result->{gin_create_time}, $gin_result->{highfreq_excluded};
    
    # DNA4 Test 4: Search with =% operator WITHOUT GIN index
    print "\nDNA4 Test 4: Search with =% operator WITHOUT GIN index (sequential scan)\n";
    my ($seq_search_time, $seq_search_result) = measure_time("dna4_seq_search_scoring", sub {
        my $table_name = "simd_test_dna4_seq_" . $simd->{value};
        $table_name =~ s/\./_/g;
        
        create_base_table($table_name, 'DNA4');
        
        # No GIN index creation - will force sequential scan
        
        # Store results for comparison
        my @all_results;
        my $batch_size = 10;  # Smaller batch size for sequential scan
        my $test_queries = 100;  # Fewer queries for sequential scan test
        
        # Test with subset of query sequences
        for (my $batch_start = 0; $batch_start < $test_queries && $batch_start < @dna4_query_sequences; $batch_start += $batch_size) {
            my $batch_end = $batch_start + $batch_size - 1;
            $batch_end = $test_queries - 1 if $batch_end >= $test_queries;
            $batch_end = $#dna4_query_sequences if $batch_end > $#dna4_query_sequences;
            
            for (my $i = $batch_start; $i <= $batch_end; $i++) {
                my $query = $dna4_query_sequences[$i];
                
                # Get matches with scores (sequential scan)
                my $sth = $dbh->prepare("
                    SELECT id, 
                           kmersearch_rawscore(seq, ?) as rawscore,
                           kmersearch_correctedscore(seq, ?) as correctedscore
                    FROM $table_name
                    WHERE seq =% ?
                    ORDER BY id
                ");
                $sth->execute($query, $query, $query);
                
                my @matches;
                while (my $row = $sth->fetchrow_hashref()) {
                    push @matches, {
                        id => $row->{id},
                        rawscore => $row->{rawscore},
                        correctedscore => $row->{correctedscore}
                    };
                }
                $sth->finish();
                
                push @all_results, {
                    query_idx => $i,
                    match_count => scalar(@matches),
                    matches => \@matches
                };
            }
            
            printf "    Processed %d/%d queries...\r", 
                   ($batch_end + 1), $test_queries;
        }
        print "\n";
        
        # Calculate summary statistics
        my $total_matches = 0;
        my $total_rawscore = 0;
        my $total_correctedscore = 0;
        my $queries_with_matches = 0;
        
        foreach my $result (@all_results) {
            $total_matches += $result->{match_count};
            $queries_with_matches++ if $result->{match_count} > 0;
            
            foreach my $match (@{$result->{matches}}) {
                $total_rawscore += $match->{rawscore};
                $total_correctedscore += $match->{correctedscore};
            }
        }
        
        # Clean up
        $dbh->do("DROP TABLE $table_name");
        
        return {
            query_count => scalar(@all_results),
            queries_with_matches => $queries_with_matches,
            total_matches => $total_matches,
            avg_matches_per_query => sprintf("%.2f", $total_matches / scalar(@all_results)),
            total_rawscore => $total_rawscore,
            total_correctedscore => $total_correctedscore,
            results_sample => [@all_results[0..9]]  # Keep first 10 for comparison
        };
    });
    
    $dna4_timing_results{$simd->{value}}{seq_search_scoring} = $seq_search_time;
    $dna4_test_results{$simd->{value}}{seq_search_scoring} = $seq_search_result;
    printf "  Time: %.3fs, Total matches: %d, Queries with matches: %d/%d\n", 
           $seq_search_time, $seq_search_result->{total_matches}, 
           $seq_search_result->{queries_with_matches}, $seq_search_result->{query_count};
    printf "  Avg matches/query: %s, Total rawscore: %d, Total correctedscore: %d\n",
           $seq_search_result->{avg_matches_per_query}, 
           $seq_search_result->{total_rawscore}, $seq_search_result->{total_correctedscore};
    
    # DNA4 Test 5: Search with =% operator WITH GIN index
    print "\nDNA4 Test 5: Search with =% operator WITH GIN index\n";
    my ($search_time, $search_result) = measure_time("dna4_search_scoring", sub {
        my $table_name = $gin_result->{table_name};
        
        # Store results for comparison
        my @all_results;
        my $batch_size = 100;  # Process queries in batches
        
        # Test with all 1000 query sequences
        for (my $batch_start = 0; $batch_start < @dna4_query_sequences; $batch_start += $batch_size) {
            my $batch_end = $batch_start + $batch_size - 1;
            $batch_end = $#dna4_query_sequences if $batch_end > $#dna4_query_sequences;
            
            for (my $i = $batch_start; $i <= $batch_end; $i++) {
                my $query = $dna4_query_sequences[$i];
                
                # Get matches with scores
                my $sth = $dbh->prepare("
                    SELECT id, 
                           kmersearch_rawscore(seq, ?) as rawscore,
                           kmersearch_correctedscore(seq, ?) as correctedscore
                    FROM $table_name
                    WHERE seq =% ?
                    ORDER BY id
                ");
                $sth->execute($query, $query, $query);
                
                my @matches;
                while (my $row = $sth->fetchrow_hashref()) {
                    push @matches, {
                        id => $row->{id},
                        rawscore => $row->{rawscore},
                        correctedscore => $row->{correctedscore}
                    };
                }
                $sth->finish();
                
                push @all_results, {
                    query_idx => $i,
                    match_count => scalar(@matches),
                    matches => \@matches
                };
            }
            
            printf "    Processed %d/%d queries...\r", 
                   ($batch_end + 1), scalar(@dna4_query_sequences);
        }
        print "\n";
        
        # Calculate summary statistics
        my $total_matches = 0;
        my $total_rawscore = 0;
        my $total_correctedscore = 0;
        my $queries_with_matches = 0;
        
        foreach my $result (@all_results) {
            $total_matches += $result->{match_count};
            $queries_with_matches++ if $result->{match_count} > 0;
            
            foreach my $match (@{$result->{matches}}) {
                $total_rawscore += $match->{rawscore};
                $total_correctedscore += $match->{correctedscore};
            }
        }
        
        # Clean up
        $dbh->do("DROP TABLE $table_name CASCADE");
        
        return {
            query_count => scalar(@all_results),
            queries_with_matches => $queries_with_matches,
            total_matches => $total_matches,
            avg_matches_per_query => sprintf("%.2f", $total_matches / scalar(@all_results)),
            total_rawscore => $total_rawscore,
            total_correctedscore => $total_correctedscore,
            results_sample => [@all_results[0..9]]  # Keep first 10 for comparison
        };
    });
    
    $dna4_timing_results{$simd->{value}}{search_scoring} = $search_time;
    $dna4_test_results{$simd->{value}}{search_scoring} = $search_result;
    printf "  Time: %.3fs, Total matches: %d, Queries with matches: %d/%d\n", 
           $search_time, $search_result->{total_matches}, 
           $search_result->{queries_with_matches}, $search_result->{query_count};
    printf "  Avg matches/query: %s, Total rawscore: %d, Total correctedscore: %d\n\n",
           $search_result->{avg_matches_per_query}, 
           $search_result->{total_rawscore}, $search_result->{total_correctedscore};
}

# Compare DNA4 results across SIMD capabilities
print "=" x 70 . "\n";
print "DNA4 Result Comparison\n";
print "=" x 70 . "\n\n";

# Compare DNA4 frequency analysis results
print "DNA4 Frequency Analysis Comparison:\n";
my $dna4_base_freq = $dna4_test_results{$simd_capabilities[0]->{value}}{frequency_analysis};
my $dna4_freq_consistent = 1;

foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
    my $freq = $dna4_test_results{$simd->{value}}{frequency_analysis};
    
    if ($freq->{highfreq_kmers_count} != $dna4_base_freq->{highfreq_kmers_count}) {
        printf "  WARNING: Different DNA4 high-freq k-mer count for %s: %d vs %d\n",
               $simd->{name}, $freq->{highfreq_kmers_count}, $dna4_base_freq->{highfreq_kmers_count};
        $dna4_freq_consistent = 0;
    }
    
    # Compare sample k-mers
    my $base_kmers = $dna4_base_freq->{highfreq_kmers_sample};
    my $curr_kmers = $freq->{highfreq_kmers_sample};
    
    my $sample_match = 1;
    for (my $i = 0; $i < @$base_kmers && $i < @$curr_kmers; $i++) {
        if ($base_kmers->[$i]{kmer} ne $curr_kmers->[$i]{kmer} ||
            $base_kmers->[$i]{reason} ne $curr_kmers->[$i]{reason}) {
            $sample_match = 0;
            last;
        }
    }
    
    if (!$sample_match) {
        print "  WARNING: DNA4 K-mer sample mismatch for $simd->{name}\n";
        $dna4_freq_consistent = 0;
    }
}

if ($dna4_freq_consistent) {
    print "  ✓ DNA4 High-frequency k-mer analysis results are consistent across all SIMD capabilities\n";
}

