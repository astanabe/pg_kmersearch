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

# Test parameters
my $sequence_length = 4096;
my $total_sequences = 10000;
my $num_query_sequences = 1000;
my $max_appearance_rate = 0.1;
my $kmer_size = 16;

# SIMD capabilities for AMD64 - using numeric values for force_simd_capability
my @simd_capabilities = (
    { name => 'None',         value => 0 },
    { name => 'AVX2',         value => 1 },
    { name => 'BMI2',         value => 2 },
    { name => 'AVX512F',      value => 3 },
    { name => 'AVX512BW',     value => 4 },
    { name => 'AVX512VBMI',   value => 5 },
    { name => 'AVX512VBMI2',  value => 6 },
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

# Select query sequences (first 1000)
my @query_sequences = @sequences[0..($num_query_sequences-1)];

print "Generated " . scalar(@sequences) . " sequences\n";
print "Selected " . scalar(@query_sequences) . " query sequences\n\n";

# Get initial SIMD capability
my ($auto_capability) = $dbh->selectrow_array("SELECT kmersearch_simd_capability()");
print "Auto-detected SIMD capability: $auto_capability\n\n";

# Store results for comparison
my %test_results;
my %timing_results;

# Function to measure execution time
sub measure_time {
    my ($name, $code) = @_;
    my $start = time();
    my $result = $code->();
    my $end = time();
    return ($end - $start, $result);
}

# Create base table function
sub create_base_table {
    my ($table_name) = @_;
    
    # Drop existing table
    $dbh->do("DROP TABLE IF EXISTS $table_name CASCADE");
    
    # Create table
    $dbh->do("CREATE TABLE $table_name (id SERIAL PRIMARY KEY, seq DNA2)");
    
    # Insert sequences
    my $sth = $dbh->prepare("INSERT INTO $table_name (seq) VALUES (?)");
    foreach my $seq (@sequences) {
        $sth->execute($seq);
    }
    $sth->finish();
}

# Additional Input/Output Tests
print "=" x 70 . "\n";
print "Additional Input/Output Tests\n";
print "=" x 70 . "\n\n";

# Test 1: Scalar input/Scalar output
print "Test 1: Scalar input (force_simd_capability=0) / Scalar output (force_simd_capability=0)\n";
{
    my $table_name = "simd_test_scalar_io";
    
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
    
    printf "  Result: %d/%d sequences match (%.1f%%)\n", 
           $match_count, $row_count, ($match_count / $row_count) * 100;
    
    # Clean up
    $dbh->do("DROP TABLE $table_name");
}

# Test 2: Scalar input/SIMD output (all SIMD values)
print "\nTest 2: Scalar input (force_simd_capability=0) / SIMD output (all force_simd_capability>0)\n";
{
    my $table_name = "simd_test_scalar_input";
    
    # Set to scalar mode for input
    $dbh->do("SET kmersearch.force_simd_capability = 0");
    
    # Create table and insert data
    create_base_table($table_name);
    
    # Test with each SIMD capability for output
    foreach my $simd (@simd_capabilities) {
        next if $simd->{value} == 0;  # Skip scalar
        
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
        
        printf "  $simd->{name} output: %d/%d sequences match (%.1f%%)\n", 
               $match_count, $row_count, ($match_count / $row_count) * 100;
    }
    
    # Clean up
    $dbh->do("DROP TABLE $table_name");
}

# Test 3: SIMD input/Scalar output (all SIMD values)
print "\nTest 3: SIMD input (all force_simd_capability>0) / Scalar output (force_simd_capability=0)\n";
{
    foreach my $simd (@simd_capabilities) {
        next if $simd->{value} == 0;  # Skip scalar
        
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
        
        printf "  $simd->{name} input: %d/%d sequences match (%.1f%%)\n", 
               $match_count, $row_count, ($match_count / $row_count) * 100;
        
        # Clean up
        $dbh->do("DROP TABLE $table_name");
    }
}

print "\n";

# Test each SIMD capability
foreach my $simd (@simd_capabilities) {
    print "=" x 70 . "\n";
    print "Testing with kmersearch.force_simd_capability = '$simd->{value}'\n";
    print "=" x 70 . "\n";
    
    # Set SIMD capability
    $dbh->do("SET kmersearch.force_simd_capability = $simd->{value}");
    $dbh->do("SET kmersearch.kmer_size = $kmer_size");
    
    # Verify setting
    my ($current) = $dbh->selectrow_array("SELECT kmersearch_simd_capability()");
    print "Current SIMD capability: $current\n\n";
    
    # Test 1: Input/Output verification
    print "Test 1: Input/Output verification\n";
    my ($io_time, $io_result) = measure_time("input_output", sub {
        my $table_name = "simd_test_io_" . $simd->{value};
        $table_name =~ s/\./_/g;  # Replace dots with underscores
        
        create_base_table($table_name);
        
        # Verify data integrity
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
    
    $timing_results{$simd->{value}}{input_output} = $io_time;
    $test_results{$simd->{value}}{input_output} = $io_result;
    printf "  Time: %.3fs, Matched: %d/%d (%.1f%%)\n", 
           $io_time, $io_result->{matched}, $io_result->{total_rows},
           ($io_result->{matched} / $io_result->{total_rows}) * 100;
    
    # Test 2: K-mer frequency analysis
    print "\nTest 2: K-mer frequency analysis (max_appearance_rate=$max_appearance_rate)\n";
    my ($freq_time, $freq_result) = measure_time("frequency_analysis", sub {
        my $table_name = "simd_test_freq_" . $simd->{value};
        $table_name =~ s/\./_/g;
        
        create_base_table($table_name);
        
        # Set parameters
        $dbh->do("SET kmersearch.max_appearance_rate = $max_appearance_rate");
        
        # Perform analysis
        my $sth = $dbh->prepare("SELECT * FROM kmersearch_perform_highfreq_analysis(?, ?)");
        $sth->execute($table_name, 'seq');
        my $analysis = $sth->fetchrow_hashref();
        $sth->finish();
        
        # Get high-frequency k-mers for comparison
        $sth = $dbh->prepare("
            SELECT kmer_encoded::text, appearance_count 
            FROM kmersearch_highfreq_kmer 
            WHERE table_oid = ?::regclass::oid 
            ORDER BY kmer_encoded::text
            LIMIT 20
        ");
        $sth->execute($table_name);
        
        my @highfreq_kmers;
        while (my ($kmer, $count) = $sth->fetchrow_array()) {
            push @highfreq_kmers, { kmer => $kmer, count => $count };
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
    
    $timing_results{$simd->{value}}{frequency_analysis} = $freq_time;
    $test_results{$simd->{value}}{frequency_analysis} = $freq_result;
    printf "  Time: %.3fs, High-freq k-mers: %d, Workers: %d\n", 
           $freq_time, $freq_result->{highfreq_kmers_count}, $freq_result->{parallel_workers};
    
    # Test 3: GIN index construction with high-frequency k-mer exclusion
    print "\nTest 3: GIN index construction with high-frequency k-mer exclusion\n";
    my ($gin_time, $gin_result) = measure_time("gin_index", sub {
        my $table_name = "simd_test_gin_" . $simd->{value};
        $table_name =~ s/\./_/g;
        
        create_base_table($table_name);
        
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
        my $gin_start = time();
        $dbh->do("CREATE INDEX $index_name ON $table_name USING gin(seq gin_dna2_kmer_ops)");
        my $gin_create_time = time() - $gin_start;
        
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
    
    $timing_results{$simd->{value}}{gin_index} = $gin_time;
    $test_results{$simd->{value}}{gin_index} = $gin_result;
    printf "  Time: %.3fs (index creation: %.3fs), High-freq k-mers excluded: %d\n", 
           $gin_time, $gin_result->{gin_create_time}, $gin_result->{highfreq_excluded};
    
    # Test 4: Search with =% operator and score calculation
    print "\nTest 4: Search with =% operator and score calculation\n";
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
                           kmersearch_rawscore(seq, ?::DNA2) as rawscore,
                           kmersearch_correctedscore(seq, ?::DNA2) as correctedscore
                    FROM $table_name
                    WHERE seq =% ?::DNA2
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

# Compare results across SIMD capabilities
print "=" x 70 . "\n";
print "Result Comparison\n";
print "=" x 70 . "\n\n";

# Compare frequency analysis results
print "Frequency Analysis Comparison:\n";
my $base_freq = $test_results{$simd_capabilities[0]->{value}}{frequency_analysis};
my $freq_consistent = 1;

foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
    my $freq = $test_results{$simd->{value}}{frequency_analysis};
    
    if ($freq->{highfreq_kmers_count} != $base_freq->{highfreq_kmers_count}) {
        printf "  WARNING: Different high-freq k-mer count for %s: %d vs %d\n",
               $simd->{name}, $freq->{highfreq_kmers_count}, $base_freq->{highfreq_kmers_count};
        $freq_consistent = 0;
    }
    
    # Compare sample k-mers
    my $base_kmers = $base_freq->{highfreq_kmers_sample};
    my $curr_kmers = $freq->{highfreq_kmers_sample};
    
    my $sample_match = 1;
    for (my $i = 0; $i < @$base_kmers && $i < @$curr_kmers; $i++) {
        if ($base_kmers->[$i]{kmer} ne $curr_kmers->[$i]{kmer} ||
            $base_kmers->[$i]{count} != $curr_kmers->[$i]{count}) {
            $sample_match = 0;
            last;
        }
    }
    
    if (!$sample_match) {
        print "  WARNING: K-mer sample mismatch for $simd->{name}\n";
        $freq_consistent = 0;
    }
}

if ($freq_consistent) {
    print "  ✓ High-frequency k-mer analysis results are consistent across all SIMD capabilities\n";
}

# Compare search results
print "\nSearch Results Comparison:\n";
my $base_search = $test_results{$simd_capabilities[0]->{value}}{search_scoring};
my $search_consistent = 1;

foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
    my $search = $test_results{$simd->{value}}{search_scoring};
    
    if ($search->{total_matches} != $base_search->{total_matches} ||
        $search->{total_rawscore} != $base_search->{total_rawscore} ||
        $search->{total_correctedscore} != $base_search->{total_correctedscore}) {
        printf "  WARNING: Different search results for %s:\n", $simd->{name};
        printf "    Matches: %d vs %d\n", $search->{total_matches}, $base_search->{total_matches};
        printf "    Rawscore: %d vs %d\n", $search->{total_rawscore}, $base_search->{total_rawscore};
        printf "    Correctedscore: %d vs %d\n", $search->{total_correctedscore}, $base_search->{total_correctedscore};
        $search_consistent = 0;
    }
}

if ($search_consistent) {
    print "  ✓ Search results are consistent across all SIMD capabilities\n";
}

# Performance summary
print "\nPerformance Summary:\n";
print "=" x 90 . "\n";
printf "%-15s %12s %12s %12s %12s %12s\n", 
       "SIMD Capability", "I/O (s)", "Freq Anal(s)", "GIN Index(s)", "Search (s)", "Total (s)";
print "-" x 90 . "\n";

foreach my $simd (@simd_capabilities) {
    my $total = $timing_results{$simd->{value}}{input_output} +
                $timing_results{$simd->{value}}{frequency_analysis} +
                $timing_results{$simd->{value}}{gin_index} +
                $timing_results{$simd->{value}}{search_scoring};
    
    printf "%-15s %12.3f %12.3f %12.3f %12.3f %12.3f\n",
        $simd->{name},
        $timing_results{$simd->{value}}{input_output},
        $timing_results{$simd->{value}}{frequency_analysis},
        $timing_results{$simd->{value}}{gin_index},
        $timing_results{$simd->{value}}{search_scoring},
        $total;
}

# Calculate speedup
if (@simd_capabilities > 1) {
    print "\nSpeedup relative to 'None':\n";
    print "=" x 90 . "\n";
    printf "%-15s %12s %12s %12s %12s %12s\n",
           "SIMD Capability", "I/O", "Freq Anal", "GIN Index", "Search", "Total";
    print "-" x 90 . "\n";
    
    my $base_value = $simd_capabilities[0]->{value};
    my $base_timing = $timing_results{$base_value};
    my $base_total = $base_timing->{input_output} + $base_timing->{frequency_analysis} +
                     $base_timing->{gin_index} + $base_timing->{search_scoring};
    
    foreach my $simd (@simd_capabilities[1..$#simd_capabilities]) {
        my $curr_total = $timing_results{$simd->{value}}{input_output} +
                         $timing_results{$simd->{value}}{frequency_analysis} +
                         $timing_results{$simd->{value}}{gin_index} +
                         $timing_results{$simd->{value}}{search_scoring};
        
        printf "%-15s %12.2fx %12.2fx %12.2fx %12.2fx %12.2fx\n",
            $simd->{name},
            $base_timing->{input_output} / $timing_results{$simd->{value}}{input_output},
            $base_timing->{frequency_analysis} / $timing_results{$simd->{value}}{frequency_analysis},
            $base_timing->{gin_index} / $timing_results{$simd->{value}}{gin_index},
            $base_timing->{search_scoring} / $timing_results{$simd->{value}}{search_scoring},
            $base_total / $curr_total;
    }
}

# Reset to auto-detect
print "\nResetting to auto-detected SIMD capability...\n";
$dbh->do("SET kmersearch.force_simd_capability = -1");
my ($final) = $dbh->selectrow_array("SELECT kmersearch_simd_capability()");
print "Final SIMD capability: $final\n";

# Disconnect
$dbh->disconnect();

print "\nAll tests completed successfully!\n";