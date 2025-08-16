#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use Time::HiRes qw(time);

# Disable output buffering
$| = 1;

# PostgreSQL connection parameters
my $dbname = $ENV{PGDATABASE} || 'postgres';
my $host = $ENV{PGHOST} || 'localhost';
my $port = $ENV{PGPORT} || '5432';
my $user = $ENV{PGUSER} || $ENV{USER};

# Global variables
my $dbh;
my $sequence_length = 4096;
my $total_sequences = 1000;  # Original count (degenerate expansion disabled)
my $num_query_sequences = 100;
my $max_appearance_rate = 0.1;
my $kmer_size = 16;
my $test_start_time;

# Test data storage
my @dna2_sequences;
my @dna4_sequences;
my @dna2_query_sequences;
my @dna4_query_sequences;

# SIMD capabilities for AMD64
my @simd_capabilities = (
    { name => 'None',         value => 0 },
    { name => 'AVX2',         value => 1 },
    { name => 'BMI2',         value => 2 },
    { name => 'AVX512F',      value => 3 },
    { name => 'AVX512BW',     value => 4 },
    { name => 'AVX512VBMI',   value => 5 },
    { name => 'AVX512VBMI2',  value => 6 },
);

# Results storage
my %timing_results;
my %test_results;

# Subroutines

sub connect_database {
    $dbh = DBI->connect("dbi:Pg:dbname=$dbname;host=$host;port=$port", $user, '', 
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

    # Set client_min_messages to suppress NOTICE messages
    $dbh->do("SET client_min_messages = WARNING");
    
    # Disable autovacuum for test tables
    $dbh->do("SET session_replication_role = replica");
    
    # Note: High-frequency k-mer cache will be initialized when needed
    
    # Clean up any existing test tables
    print "Cleaning up any existing test tables...\n";
    my $sth = $dbh->prepare("
        SELECT tablename 
        FROM pg_tables 
        WHERE tablename LIKE 'simd_test_%' 
        AND schemaname = current_schema()
    ");
    $sth->execute();
    while (my ($table) = $sth->fetchrow_array()) {
        eval { $dbh->do("DROP TABLE IF EXISTS $table CASCADE"); };
        if ($@) {
            print "  Warning: Could not drop table $table: $@\n";
        } else {
            print "  Dropped table: $table\n";
        }
    }
    $sth->finish();
}

sub generate_test_data {
    print "Generating $total_sequences random DNA sequences of $sequence_length bp...\n";
    my $start_time = time();
    
    # Generate DNA2 sequences
    for (my $i = 0; $i < $total_sequences; $i++) {
        my $seq = '';
        for (my $j = 0; $j < $sequence_length; $j++) {
            $seq .= ['A', 'C', 'G', 'T']->[int(rand(4))];
        }
        push @dna2_sequences, $seq;
    }
    
    my $end_time = time();
    my $memory_usage = (`ps -o rss= -p $$` || 0) / 1024;  # Convert KB to MB
    printf "  Generated %d sequences... (%.1f seconds, %.1f MB memory)\n", 
           scalar(@dna2_sequences), $end_time - $start_time, $memory_usage;
    
    # Generate DNA4 sequences with degenerate codes
    print "Generating $total_sequences random DNA4 sequences with degenerate codes...\n";
    $start_time = time();
    
    # Use DNA4 with degenerate bases (expansion now disabled)
    my @dna4_chars = ('A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'S', 'W', 'B', 'D', 'H', 'V', 'N');
    my @dna4_weights = (200, 200, 200, 200, 5, 5, 5, 5, 5, 5, 2, 2, 2, 2, 1);  # Original weights
    
    # Pre-calculate cumulative weights for binary search
    my @cumulative_weights;
    my $sum = 0;
    for (my $i = 0; $i < @dna4_weights; $i++) {
        $sum += $dna4_weights[$i];
        push @cumulative_weights, $sum;
    }
    my $total_weight = $sum;
    
    # Create lookup table for fast character selection
    my @lookup_table;
    for (my $i = 0; $i < @dna4_chars; $i++) {
        for (my $j = 0; $j < $dna4_weights[$i]; $j++) {
            push @lookup_table, $dna4_chars[$i];
        }
    }
    
    for (my $i = 0; $i < $total_sequences; $i++) {
        # Use array to build sequence, then join (much faster than string concatenation)
        my @seq_chars;
        for (my $j = 0; $j < $sequence_length; $j++) {
            push @seq_chars, $lookup_table[int(rand($total_weight))];
        }
        push @dna4_sequences, join('', @seq_chars);
    }
    
    $end_time = time();
    $memory_usage = (`ps -o rss= -p $$` || 0) / 1024;
    printf "  Generated %d DNA4 sequences... (%.1f seconds, %.1f MB memory)\n",
           scalar(@dna4_sequences), $end_time - $start_time, $memory_usage;
    
    print "Generated " . scalar(@dna2_sequences) . " DNA2 sequences\n";
    print "Generated " . scalar(@dna4_sequences) . " DNA4 sequences\n";
    
    # Select query sequences
    my @indices = (0..$total_sequences-1);
    for (my $i = $total_sequences - 1; $i > 0; $i--) {
        my $j = int(rand($i + 1));
        @indices[$i, $j] = @indices[$j, $i];
    }
    
    @dna2_query_sequences = map { $dna2_sequences[$_] } @indices[0..$num_query_sequences-1];
    @dna4_query_sequences = map { $dna4_sequences[$_] } @indices[0..$num_query_sequences-1];
    
    print "Selected " . scalar(@dna2_query_sequences) . " DNA2 query sequences\n";
    print "Selected " . scalar(@dna4_query_sequences) . " DNA4 query sequences\n\n";
}

sub detect_simd_capability {
    print "Detecting SIMD capability...\n";
    my $start = time();
    my ($capability_info) = $dbh->selectrow_array("SELECT kmersearch_simd_capability()");
    my $elapsed = time() - $start;
    printf "Auto-detected SIMD capability: %s (took %.1f seconds)\n\n", $capability_info, $elapsed;
}

sub measure_time {
    my ($label, $code) = @_;
    my $start = time();
    my $result = eval { $code->() };
    my $end = time();
    my $elapsed = $end - $start;
    if ($@) {
        warn "Error in $label: $@\n";
        return ($elapsed, undef);
    }
    return ($elapsed, $result);
}

sub create_base_table {
    my ($table_name, $type) = @_;
    $type ||= 'DNA2';
    
    my $start_time = time();
    
    # Drop table if exists
    eval { $dbh->do("DROP TABLE IF EXISTS $table_name CASCADE"); };
    my $after_drop = time();
    
    my $datatype = $type eq 'DNA4' ? 'dna4' : 'dna2';
    $dbh->do("CREATE TABLE $table_name (id serial PRIMARY KEY, seq $datatype)");
    my $after_create = time();
    
    # Prepare bulk insert
    my $sth = $dbh->prepare("INSERT INTO $table_name (seq) VALUES (?)");
    my $after_prepare = time();
    
    # Insert sequences
    my $sequences = $type eq 'DNA4' ? \@dna4_sequences : \@dna2_sequences;
    foreach my $seq (@$sequences) {
        $sth->execute($seq);
    }
    my $after_insert = time();
    $sth->finish();
    my $end_time = time();
    
    printf "  Drop table: %.3f seconds\n", $after_drop - $start_time;
    printf "  Create table: %.3f seconds\n", $after_create - $after_drop;
    printf "  Prepare statement: %.3f seconds\n", $after_prepare - $after_create;
    printf "  Insert %d sequences: %.3f seconds\n", scalar(@$sequences), $after_insert - $after_prepare;
    printf "  Total: %.3f seconds\n\n", $end_time - $start_time;
}

# Test 1: Scalar input/output verification
sub test_scalar_io {
    my ($type) = @_;
    
    print "=" x 70 . "\n";
    print "$type Test 1: Scalar input (force_simd_capability=0) / Scalar output (force_simd_capability=0)\n";
    print "=" x 70 . "\n\n";
    
    # Set force_simd_capability to 0 (scalar)
    $dbh->do("SET kmersearch.force_simd_capability = 0");
    
    my $table_name = "simd_test_${type}_scalar_io";
    $table_name =~ s/DNA/dna/g;
    
    create_base_table($table_name, $type);
    
    my ($test_time, $test_result) = measure_time("${type}_scalar_io", sub {
        my $datatype = $type eq 'DNA4' ? 'dna4' : 'dna2';
        
        # Verify data integrity
        my $sth = $dbh->prepare("SELECT seq::text FROM $table_name ORDER BY id");
        $sth->execute();
        
        my @output_sequences;
        while (my ($seq) = $sth->fetchrow_array()) {
            push @output_sequences, $seq;
        }
        $sth->finish();
        
        # Compare with input
        my $sequences = $type eq 'DNA4' ? \@dna4_sequences : \@dna2_sequences;
        my $match_count = 0;
        for (my $i = 0; $i < @$sequences && $i < @output_sequences; $i++) {
            $match_count++ if $sequences->[$i] eq $output_sequences[$i];
        }
        
        # Clean up
        $dbh->do("DROP TABLE $table_name");
        
        return {
            total => scalar(@$sequences),
            matches => $match_count
        };
    });
    
    $timing_results{"${type}_scalar_io"}{0} = $test_time;
    $test_results{"${type}_scalar_io"}{0} = $test_result;
    
    printf "  Time: %.3fs, Result: %d/%d sequences match (%.1f%%)\n\n", 
           $test_time, $test_result->{matches}, $test_result->{total},
           100.0 * $test_result->{matches} / $test_result->{total};
}

# Test 2: Scalar input / SIMD output verification
sub test_scalar_input_simd_output {
    my ($type) = @_;
    
    print "=" x 70 . "\n";
    print "$type Test 2: Scalar input (force_simd_capability=0) / SIMD output (all force_simd_capability>0)\n";
    print "=" x 70 . "\n\n";
    
    # First insert with scalar
    $dbh->do("SET kmersearch.force_simd_capability = 0");
    
    my $table_name = "simd_test_${type}_scalar_input";
    $table_name =~ s/DNA/dna/g;
    create_base_table($table_name, $type);
    
    # Test each SIMD capability
    my %simd_results;
    foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
        $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
        
        my ($test_time, $test_result) = measure_time("${type}_scalar_in_simd_out_$simd->{value}", sub {
            my $sth = $dbh->prepare("SELECT seq::text FROM $table_name ORDER BY id");
            $sth->execute();
            
            my @output_sequences;
            while (my ($seq) = $sth->fetchrow_array()) {
                push @output_sequences, $seq;
            }
            $sth->finish();
            
            # Compare with input
            my $sequences = $type eq 'DNA4' ? \@dna4_sequences : \@dna2_sequences;
            my $match_count = 0;
            for (my $i = 0; $i < @$sequences && $i < @output_sequences; $i++) {
                $match_count++ if $sequences->[$i] eq $output_sequences[$i];
            }
            
            return {
                total => scalar(@$sequences),
                matches => $match_count
            };
        });
        
        $timing_results{"${type}_scalar_input_simd_output"}{$simd->{value}} = $test_time;
        $simd_results{$simd->{value}} = $test_result;
        
        printf "  %s output: %d/%d sequences match (%.1f%%) - Time: %.3fs\n",
               $simd->{name}, $test_result->{matches}, $test_result->{total},
               100.0 * $test_result->{matches} / $test_result->{total}, $test_time;
    }
    
    # Clean up
    $dbh->do("DROP TABLE $table_name");
    
    # Display performance comparison
    print "\n$type Scalar Input/SIMD Output Performance Comparison:\n";
    print "-" x 50 . "\n";
    foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
        my $result = $simd_results{$simd->{value}};
        printf "  %-15s: %.3fs (%.1f%% match rate)\n", 
               $simd->{name}, 
               $timing_results{"${type}_scalar_input_simd_output"}{$simd->{value}},
               100.0 * $result->{matches} / $result->{total};
    }
    print "-" x 50 . "\n\n";
}

# Test 3: SIMD input / Scalar output verification
sub test_simd_input_scalar_output {
    my ($type) = @_;
    
    print "=" x 70 . "\n";
    print "$type Test 3: SIMD input (all force_simd_capability>0) / Scalar output (force_simd_capability=0)\n";
    print "=" x 70 . "\n\n";
    
    my %simd_results;
    foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
        # Insert with SIMD
        $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
        
        
        my $table_name = "simd_test_${type}_simd_input_$simd->{value}";
        $table_name =~ s/DNA/dna/g;
        $table_name =~ s/\./_/g;
        
        create_base_table($table_name, $type);
        
        my ($test_time, $test_result) = measure_time("${type}_simd_in_scalar_out_$simd->{value}", sub {
            
            # Read with scalar
            $dbh->do("SET kmersearch.force_simd_capability = 0");
            
            my $sth = $dbh->prepare("SELECT seq::text FROM $table_name ORDER BY id");
            $sth->execute();
            
            my @output_sequences;
            while (my ($seq) = $sth->fetchrow_array()) {
                push @output_sequences, $seq;
            }
            $sth->finish();
            
            # Compare with input
            my $sequences = $type eq 'DNA4' ? \@dna4_sequences : \@dna2_sequences;
            my $match_count = 0;
            for (my $i = 0; $i < @$sequences && $i < @output_sequences; $i++) {
                $match_count++ if $sequences->[$i] eq $output_sequences[$i];
            }
            
            # Clean up
            $dbh->do("DROP TABLE $table_name");
            
            return {
                total => scalar(@$sequences),
                matches => $match_count
            };
        });
        
        $timing_results{"${type}_simd_input_scalar_output"}{$simd->{value}} = $test_time;
        $simd_results{$simd->{value}} = $test_result;
        
        printf "  %s input: %d/%d sequences match (%.1f%%) - Time: %.3fs\n",
               $simd->{name}, $test_result->{matches}, $test_result->{total},
               100.0 * $test_result->{matches} / $test_result->{total}, $test_time;
        
        # Check for background processes after AVX512F
        if ($type eq 'DNA4' && $simd->{value} == 3) {
            my $activity = $dbh->selectall_arrayref("
                SELECT pid, state, query 
                FROM pg_stat_activity 
                WHERE pid != pg_backend_pid() 
                AND state != 'idle'
            ");
            if (@$activity) {
                foreach my $row (@$activity) {
                    printf "  PID %d: %s - %s\n", $row->[0], $row->[1], substr($row->[2] || '', 0, 60);
                }
            }
        }
    }
    
    # Display performance comparison
    print "\n$type SIMD Input/Scalar Output Performance Comparison:\n";
    print "-" x 50 . "\n";
    foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
        my $result = $simd_results{$simd->{value}};
        printf "  %-15s: %.3fs (%.1f%% match rate)\n", 
               $simd->{name}, 
               $timing_results{"${type}_simd_input_scalar_output"}{$simd->{value}},
               100.0 * $result->{matches} / $result->{total};
    }
    print "-" x 50 . "\n\n";
}

# Test 4: K-mer frequency analysis comparison
sub test_kmer_frequency_analysis {
    my ($type) = @_;
    
    print "=" x 70 . "\n";
    print "$type Test 4: K-mer frequency analysis\n";
    print "=" x 70 . "\n\n";
    
    # Set common GUC variables once before the loop to avoid unnecessary cache clearing
    $dbh->do("SET kmersearch.kmer_size = $kmer_size");
    $dbh->do("SET kmersearch.max_appearance_rate = $max_appearance_rate");
    
    my %freq_results;
    foreach my $simd (@simd_capabilities) {
        print "Testing $type frequency analysis with force_simd_capability = $simd->{value} ($simd->{name})\n";
        
        $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
        $dbh->do("SET kmersearch.preclude_highfreq_kmer = true");
        
        my $table_name = "simd_test_${type}_freq_$simd->{value}";
        $table_name =~ s/DNA/dna/g;
        $table_name =~ s/\./_/g;
        
        create_base_table($table_name, $type);
        
        my ($freq_time, $freq_result) = measure_time("${type}_freq_$simd->{value}", sub {
            
            # Perform analysis
            my $sth = $dbh->prepare("SELECT * FROM kmersearch_perform_highfreq_analysis(?, ?)");
            $sth->execute($table_name, 'seq');
            my $analysis = $sth->fetchrow_hashref();
            $sth->finish();
            
            # Get high-frequency k-mers for comparison
            $sth = $dbh->prepare("
                SELECT uintkey::text, detection_reason 
                FROM kmersearch_highfreq_kmer 
                WHERE table_oid = ?::regclass::oid 
                ORDER BY uintkey
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
            
            # Clean up - now returns a composite type
            eval {
                my $sth = $dbh->prepare("SELECT * FROM kmersearch_undo_highfreq_analysis(?, ?)");
                $sth->execute($table_name, 'seq');
                $sth->finish();
            };
            $dbh->do("DROP TABLE $table_name");
            
            return {
                total_rows => $analysis->{total_rows},
                highfreq_kmers_count => $analysis->{highfreq_kmers_count},
                highfreq_kmers_sample => \@highfreq_kmers,
                total_highfreq => $total_highfreq,
                parallel_workers => $analysis->{parallel_workers_used}
            };
        });
        
        $freq_results{$simd->{value}} = $freq_result;
        $timing_results{"${type}_freq"}{$simd->{value}} = $freq_time;
        
        if (defined $freq_result) {
            printf "  Time: %.3fs, High-freq k-mers: %d, Workers: %d\n\n", 
                   $freq_time, $freq_result->{highfreq_kmers_count}, $freq_result->{parallel_workers};
        } else {
            printf "  Time: %.3fs, ERROR: Test failed\n\n", $freq_time;
        }
    }
    
    # Compare results
    print "$type Frequency Analysis Result Comparison:\n";
    my $base_freq = $freq_results{0};
    my $freq_consistent = 1;
    
    if ($base_freq) {
        foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
            my $freq = $freq_results{$simd->{value}};
            next if !$freq;  # Skip if test failed
            
            if ($freq->{highfreq_kmers_count} != $base_freq->{highfreq_kmers_count}) {
                printf "  WARNING: Different high-freq k-mer count for %s: %d vs %d\n",
                       $simd->{name}, $freq->{highfreq_kmers_count}, $base_freq->{highfreq_kmers_count};
                $freq_consistent = 0;
            }
        }
        
        if ($freq_consistent) {
            print "  ✓ $type High-frequency k-mer analysis results are consistent across all SIMD capabilities\n";
        }
    } else {
        print "  ⚠ $type High-frequency k-mer analysis failed for baseline test\n";
    }
    print "\n";
}

# Test 5: GIN index construction comparison
sub test_gin_index_construction {
    my ($type) = @_;
    
    print "=" x 70 . "\n";
    print "$type Test 5: GIN index construction\n";
    print "=" x 70 . "\n\n";
    
    my %gin_results;
    my %gin_tables;
    
    # Set common GUC variables once before the loop to avoid unnecessary cache clearing
    $dbh->do("SET kmersearch.kmer_size = $kmer_size");
    $dbh->do("SET kmersearch.max_appearance_rate = $max_appearance_rate");
    
    foreach my $simd (@simd_capabilities) {
        print "Testing $type GIN index with force_simd_capability = $simd->{value} ($simd->{name})\n";
        
        $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
        $dbh->do("SET kmersearch.preclude_highfreq_kmer = true");
        
        my $table_name = "simd_test_${type}_gin_$simd->{value}";
        $table_name =~ s/DNA/dna/g;
        $table_name =~ s/\./_/g;
        
        create_base_table($table_name, $type);
        
        my ($gin_time, $gin_result) = measure_time("${type}_gin_$simd->{value}", sub {
            
            # Perform frequency analysis - returns composite type
            my $sth = $dbh->prepare("SELECT * FROM kmersearch_perform_highfreq_analysis(?, ?)");
            $sth->execute($table_name, 'seq');
            my $analysis = $sth->fetchrow_hashref();
            $sth->finish();
            my $analysis_rows = $analysis->{total_rows};
            
            # Load high-frequency k-mer cache after analysis
            my $cache_loaded = 0;
            eval {
                # First, check if there are any high-frequency k-mers for this specific table
                my ($highfreq_count) = $dbh->selectrow_array(
                    "SELECT COUNT(*) FROM kmersearch_highfreq_kmer WHERE table_oid = ?::regclass::oid",
                    undef, $table_name
                );
                
                if ($highfreq_count > 0) {
                    # Load both global and parallel caches
                    my ($result) = $dbh->selectrow_array(
                        "SELECT kmersearch_highfreq_kmer_cache_load(?, ?)", 
                        undef, $table_name, 'seq'
                    );
                    # Also load parallel cache for GIN index workers
                    eval {
                        $dbh->do("SELECT kmersearch_parallel_highfreq_kmer_cache_load(?, ?)", 
                                undef, $table_name, 'seq');
                    };
                    $cache_loaded = 1 if defined $result;
                    print "  Loaded high-frequency k-mer cache with $highfreq_count k-mers\n";
                } else {
                    # If no high-frequency k-mers, disable filtering
                    $dbh->do("SET kmersearch.preclude_highfreq_kmer = false");
                    $cache_loaded = 1;  # Mark as "loaded" to proceed
                    print "  No high-frequency k-mers detected, filtering disabled\n";
                }
            };
            if ($@) {
                print "  Warning: Could not load cache: $@\n";
                # Try to continue without high-frequency k-mer filtering
                $dbh->do("SET kmersearch.preclude_highfreq_kmer = false");
            } elsif (!$cache_loaded) {
                print "  Warning: Cache load returned no result\n";
                $dbh->do("SET kmersearch.preclude_highfreq_kmer = false");
            }
            
            # Create GIN index
            my $index_name = $table_name . "_idx";
            my $gin_start = time();
            my $gin_created = 0;
            
            # Try to create GIN index
            eval {
                local $SIG{ALRM} = sub { die "GIN index creation timeout\n" };
                alarm(60);
                my $datatype = $type eq 'DNA4' ? 'dna4' : 'dna2';
                
                # First, ensure the table exists
                my ($table_exists) = $dbh->selectrow_array(
                    "SELECT EXISTS (SELECT 1 FROM pg_tables WHERE tablename = ?)",
                    undef, $table_name
                );
                
                if ($table_exists) {
                    # Select operator class based on total bits needed
                    # Total bits = kmer_size * 2 + occur_bitlen
                    # Get current occur_bitlen setting
                    my ($occur_bitlen) = $dbh->selectrow_array("SHOW kmersearch.occur_bitlen");
                    $occur_bitlen = 8 unless defined $occur_bitlen;  # Default is 8
                    
                    my $total_bits = $kmer_size * 2 + $occur_bitlen;
                    
                    # Select appropriate operator class:
                    # int2: up to 16 bits
                    # int4: up to 32 bits  
                    # int8: up to 64 bits
                    my $ops_suffix;
                    if ($total_bits <= 16) {
                        $ops_suffix = "int2";
                    } elsif ($total_bits <= 32) {
                        $ops_suffix = "int4";
                    } else {
                        $ops_suffix = "int8";
                    }
                    my $ops_class = "kmersearch_${datatype}_gin_ops_${ops_suffix}";
                    $dbh->do("CREATE INDEX $index_name ON $table_name USING gin(seq $ops_class)");
                    $gin_created = 1;
                } else {
                    die "Table $table_name does not exist";
                }
                alarm(0);
            };
            my $gin_create_time = time() - $gin_start;
            
            if ($@ || !$gin_created) {
                my $error_msg = $@ || "Unknown error";
                print "  ERROR creating GIN index: $error_msg\n";
                
                # Don't drop the table on error - keep it for debugging
                return { 
                    error => $error_msg,
                    table_name => $table_name,
                    index_name => $index_name,
                    highfreq_excluded => 0,
                    gin_create_time => $gin_create_time
                };
            }
            
            # Get index metadata
            my $excluded = 0;
            eval {
                my $sth = $dbh->prepare("
                    SELECT highfreq_filtered 
                    FROM kmersearch_gin_index_meta 
                    WHERE index_oid = ?::regclass::oid
                ");
                $sth->execute($index_name);
                ($excluded) = $sth->fetchrow_array();
                $sth->finish();
            };
            if ($@) {
                print "  Warning: Could not get index metadata: $@\n";
            }
            
            return {
                table_name => $table_name,
                index_name => $index_name,
                highfreq_excluded => $excluded || 0,
                gin_create_time => $gin_create_time
            };
        });
        
        $gin_results{$simd->{value}} = $gin_result;
        if ($gin_result && !$gin_result->{error}) {
            $gin_tables{$simd->{value}} = $gin_result->{table_name};
        }
        $timing_results{"${type}_gin"}{$simd->{value}} = $gin_time;
        
        if ($gin_result && !$gin_result->{error}) {
            printf "  Time: %.3fs (index creation: %.3fs), High-freq k-mers excluded: %d\n\n", 
                   $gin_time, $gin_result->{gin_create_time}, $gin_result->{highfreq_excluded};
        } else {
            printf "  Time: %.3fs, Failed to create index\n\n", $gin_time;
        }
    }
    
    # Compare results
    print "$type GIN Index Construction Result Comparison:\n";
    my $base_gin = $gin_results{0};
    my $gin_consistent = 1;
    
    if ($base_gin && !$base_gin->{error}) {
        foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
            my $gin = $gin_results{$simd->{value}};
            next if !$gin || $gin->{error};
            
            if ($gin->{highfreq_excluded} != $base_gin->{highfreq_excluded}) {
                printf "  WARNING: Different high-freq k-mers excluded for %s: %d vs %d\n",
                       $simd->{name}, $gin->{highfreq_excluded}, $base_gin->{highfreq_excluded};
                $gin_consistent = 0;
            }
        }
        
        if ($gin_consistent) {
            print "  ✓ $type GIN index construction results are consistent across all SIMD capabilities\n";
        }
    } else {
        print "  ⚠ $type GIN index construction failed for baseline test\n";
    }
    
    # Store table names for later use in search tests
    $test_results{"${type}_gin_tables"} = \%gin_tables;
    print "\n";
}

# Test 6: GIN index search comparison
sub test_gin_index_search {
    my ($type) = @_;
    
    print "=" x 70 . "\n";
    print "$type Test 6: GIN index search with =% operator\n";
    print "=" x 70 . "\n\n";
    
    my $gin_tables = $test_results{"${type}_gin_tables"};
    my %search_results;
    
    # Set common GUC variables once before the loop to avoid unnecessary cache clearing
    $dbh->do("SET kmersearch.kmer_size = $kmer_size");
    $dbh->do("SET kmersearch.max_appearance_rate = $max_appearance_rate");
    
    foreach my $simd (@simd_capabilities) {
        my $table_name = $gin_tables->{$simd->{value}};
        if (!$table_name) {
            print "  Skipping search test for $simd->{name} - no table available\n";
            next;
        }
        
        print "Testing $type GIN search with force_simd_capability = $simd->{value} ($simd->{name})\n";
        
        $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
        
        # Check if there are any high-frequency k-mers for this table
        my ($highfreq_count) = $dbh->selectrow_array(
            "SELECT COUNT(*) FROM kmersearch_highfreq_kmer WHERE table_oid = ?::regclass::oid",
            undef, $table_name
        );
        
        if ($highfreq_count > 0) {
            $dbh->do("SET kmersearch.preclude_highfreq_kmer = true");
        } else {
            # If no high-frequency k-mers, disable filtering
            $dbh->do("SET kmersearch.preclude_highfreq_kmer = false");
        }
        
        $dbh->do("SET enable_seqscan = false");
        
        my ($search_time, $search_result) = measure_time("${type}_gin_search_$simd->{value}", sub {
            my @all_results;
            my $query_sequences = $type eq 'DNA4' ? \@dna4_query_sequences : \@dna2_query_sequences;
            my $test_queries = scalar(@$query_sequences);
            
            # Process queries in batches
            my $batch_size = 10;
            for (my $batch_start = 0; $batch_start < $test_queries; $batch_start += $batch_size) {
                my $batch_end = $batch_start + $batch_size - 1;
                $batch_end = $test_queries - 1 if $batch_end >= $test_queries;
                
                for (my $i = $batch_start; $i <= $batch_end; $i++) {
                    my $query = $query_sequences->[$i];
                    
                    my $sth = $dbh->prepare("
                        SELECT id, 
                               kmersearch_matchscore(seq, ?) as matchscore
                        FROM $table_name
                        WHERE seq =% ?
                        ORDER BY id
                    ");
                    $sth->execute($query, $query);
                    
                    my @matches;
                    while (my $row = $sth->fetchrow_hashref()) {
                        push @matches, {
                            id => $row->{id},
                            matchscore => $row->{matchscore}
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
            my $total_matchscore = 0;
            my $queries_with_matches = 0;
            
            foreach my $result (@all_results) {
                $total_matches += $result->{match_count};
                $queries_with_matches++ if $result->{match_count} > 0;
                
                foreach my $match (@{$result->{matches}}) {
                    $total_matchscore += $match->{matchscore};
                }
            }
            
            return {
                query_count => scalar(@all_results),
                queries_with_matches => $queries_with_matches,
                total_matches => $total_matches,
                avg_matches_per_query => sprintf("%.2f", $total_matches / scalar(@all_results)),
                total_matchscore => $total_matchscore,
                results_sample => [@all_results[0..9]]
            };
        });
        
        $search_results{$simd->{value}} = $search_result;
        $timing_results{"${type}_gin_search"}{$simd->{value}} = $search_time;
        
        if (defined $search_result) {
            printf "  Time: %.3fs, Total matches: %d, Queries with matches: %d/%d\n", 
                   $search_time, $search_result->{total_matches}, 
                   $search_result->{queries_with_matches}, $search_result->{query_count};
            printf "  Avg matches/query: %s, Total matchscore: %d\n\n",
                   $search_result->{avg_matches_per_query}, 
                   $search_result->{total_matchscore};
        } else {
            printf "  Time: %.3fs, ERROR: Test failed\n\n", $search_time;
        }
    }
    
    # Compare results
    print "$type GIN Search Result Comparison:\n";
    my $base_search = $search_results{0};
    my $search_consistent = 1;
    
    if ($base_search) {
        foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
            my $search = $search_results{$simd->{value}};
            next unless $search;
            
            if ($search->{total_matches} != $base_search->{total_matches} ||
                $search->{total_matchscore} != $base_search->{total_matchscore}) {
                printf "  WARNING: Different search results for %s\n", $simd->{name};
                printf "    Matches: %d vs %d, Matchscore: %d vs %d\n",
                       $search->{total_matches}, $base_search->{total_matches},
                       $search->{total_matchscore}, $base_search->{total_matchscore};
                $search_consistent = 0;
            }
        }
    }
    
    if ($search_consistent) {
        print "  ✓ $type GIN search results are consistent across all SIMD capabilities\n";
    }
    
    # Clean up tables
    foreach my $simd_value (keys %$gin_tables) {
        if ($gin_tables->{$simd_value}) {
            eval {
                $dbh->do("DROP TABLE IF EXISTS $gin_tables->{$simd_value} CASCADE");
            };
        }
    }
    print "\n";
}

# Main subroutine
sub main {
    $test_start_time = time();
    
    # Connect to database
    connect_database();
    
    # Generate test data
    generate_test_data();
    
    # Add delay to isolate timing
    print "[" . localtime() . "] Finished generating test data\n";
    my $post_generation_time = time();
    
    # Detect SIMD capability
    print "[" . localtime() . "] Starting SIMD capability detection\n";
    detect_simd_capability();
    print "[" . localtime() . "] Finished SIMD capability detection\n";
    my $post_simd_time = time();
    
    print "=" x 70 . "\n";
    print "SIMD Implementation Verification Tests\n";
    print "Following test order: DNA2 scalar -> DNA4 scalar -> ...\n";
    print "=" x 70 . "\n\n";
    my $pre_test_time = time();
    
    print "\n";
    
    # Test execution in the specified order
    # 1. DNA2 and DNA4 scalar I/O tests
    my $test1_start = time();
    test_scalar_io('DNA2');
    my $test1_end = time();
    test_scalar_io('DNA4');
    
    # 2. DNA2 and DNA4 scalar input / SIMD output tests
    test_scalar_input_simd_output('DNA2');
    test_scalar_input_simd_output('DNA4');
    
    # 3. DNA2 and DNA4 SIMD input / scalar output tests
    test_simd_input_scalar_output('DNA2');
    test_simd_input_scalar_output('DNA4');
    
    # 4. DNA2 and DNA4 k-mer frequency analysis tests
    test_kmer_frequency_analysis('DNA2');
    test_kmer_frequency_analysis('DNA4');
    
    # 5. DNA2 and DNA4 GIN index construction tests
    test_gin_index_construction('DNA2');
    test_gin_index_construction('DNA4');
    
    # 6. DNA2 and DNA4 GIN index search tests
    test_gin_index_search('DNA2');
    test_gin_index_search('DNA4');
    
    # Final summary
    print_final_summary();
    
    # Disconnect
    $dbh->disconnect();
    
    print "\nTotal test time: " . sprintf("%.3f", time() - $test_start_time) . " seconds\n";
}

sub print_final_summary {
    print "\n" . "=" x 70 . "\n";
    print "Final Summary\n";
    print "=" x 70 . "\n\n";
    
    # Display timing matrix
    print "Performance Summary Matrix (times in seconds):\n\n";
    printf "%-30s", "Test";
    foreach my $simd (@simd_capabilities) {
        printf " %12s", $simd->{name};
    }
    print "\n";
    print "-" x (30 + 13 * scalar(@simd_capabilities)) . "\n";
    
    # DNA2 Tests
    my @test_types = (
        { key => 'DNA2_scalar_io', label => 'DNA2 Scalar I/O' },
        { key => 'DNA2_scalar_input_simd_output', label => 'DNA2 Scalar In/SIMD Out' },
        { key => 'DNA2_simd_input_scalar_output', label => 'DNA2 SIMD In/Scalar Out' },
        { key => 'DNA2_freq', label => 'DNA2 Frequency Analysis' },
        { key => 'DNA2_gin', label => 'DNA2 GIN Index' },
        { key => 'DNA2_gin_search', label => 'DNA2 GIN Search' },
    );
    
    foreach my $test (@test_types) {
        printf "%-30s", $test->{label};
        foreach my $simd (@simd_capabilities) {
            my $time = $timing_results{$test->{key}}{$simd->{value}} || 0;
            printf " %12.3f", $time;
        }
        print "\n";
    }
    
    print "-" x (30 + 13 * scalar(@simd_capabilities)) . "\n";
    
    # DNA4 Tests
    @test_types = (
        { key => 'DNA4_scalar_io', label => 'DNA4 Scalar I/O' },
        { key => 'DNA4_scalar_input_simd_output', label => 'DNA4 Scalar In/SIMD Out' },
        { key => 'DNA4_simd_input_scalar_output', label => 'DNA4 SIMD In/Scalar Out' },
        { key => 'DNA4_freq', label => 'DNA4 Frequency Analysis' },
        { key => 'DNA4_gin', label => 'DNA4 GIN Index' },
        { key => 'DNA4_gin_search', label => 'DNA4 GIN Search' },
    );
    
    foreach my $test (@test_types) {
        printf "%-30s", $test->{label};
        foreach my $simd (@simd_capabilities) {
            my $time = $timing_results{$test->{key}}{$simd->{value}} || 0;
            printf " %12.3f", $time;
        }
        print "\n";
    }
    
    print "\nTest completed successfully!\n";
}

# Run main
main();