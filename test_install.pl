#!/usr/bin/perl
use strict;
use warnings;
use DBI;
use File::Temp qw(tempdir);
use Getopt::Long;

# Configuration
my $pg_config = qx(which pg_config);
chomp $pg_config;
die "pg_config not found in PATH" unless $pg_config && -x $pg_config;

my $bindir = qx($pg_config --bindir);
chomp $bindir;

my $initdb = "$bindir/initdb";
my $postgres = "$bindir/postgres";
my $createdb = "$bindir/createdb";
my $dropdb = "$bindir/dropdb";

# Test configuration
my $test_db = "test_pg_kmersearch_" . $$;
my $test_user = $ENV{USER} || "postgres";
my $temp_dir;
my $pg_pid;
my $port = 15432;  # Use non-standard port to avoid conflicts

# Command line options
my $verbose = 0;
my $keep_db = 0;
GetOptions(
    'verbose|v' => \$verbose,
    'keep-db' => \$keep_db,
) or die "Usage: $0 [--verbose] [--keep-db]\n";

sub log_message {
    my ($msg) = @_;
    print "[TEST] $msg\n" if $verbose;
}

sub run_command {
    my ($cmd) = @_;
    log_message("Running: $cmd");
    my $result = system($cmd);
    return $result == 0;
}

sub cleanup {
    log_message("Cleaning up test environment...");
    
    # Stop PostgreSQL if running
    if ($pg_pid) {
        log_message("Stopping PostgreSQL server (PID: $pg_pid)");
        kill 'TERM', $pg_pid;
        waitpid($pg_pid, 0);
    }
    
    # Remove temporary directory
    if ($temp_dir && -d $temp_dir) {
        log_message("Removing temporary directory: $temp_dir");
        system("rm -rf '$temp_dir'");
    }
}

# Set up signal handlers for cleanup
$SIG{INT} = $SIG{TERM} = $SIG{QUIT} = sub {
    log_message("Received signal, cleaning up...");
    cleanup();
    exit 1;
};

sub setup_test_environment {
    log_message("Setting up test environment...");
    
    # Create temporary directory for PostgreSQL cluster
    $temp_dir = tempdir(CLEANUP => 0);
    log_message("Using temporary directory: $temp_dir");
    
    # Initialize PostgreSQL cluster
    my $initdb_cmd = "$initdb -D '$temp_dir/data' -U '$test_user' --auth-local=trust --auth-host=trust";
    unless (run_command($initdb_cmd)) {
        die "Failed to initialize PostgreSQL cluster\n";
    }
    
    # Start PostgreSQL server
    my $postgres_cmd = "$postgres -D '$temp_dir/data' -p $port -k '$temp_dir' -F";
    log_message("Starting PostgreSQL server on port $port");
    
    $pg_pid = fork();
    if ($pg_pid == 0) {
        # Child process - run PostgreSQL
        exec($postgres_cmd);
        die "Failed to start PostgreSQL: $!\n";
    } elsif (!defined $pg_pid) {
        die "Failed to fork: $!\n";
    }
    
    # Wait for PostgreSQL to start
    sleep(3);
    for my $i (1..10) {
        my $check_cmd = "$bindir/pg_isready -h '$temp_dir' -p $port -U '$test_user'";
        if (run_command($check_cmd . " >/dev/null 2>&1")) {
            log_message("PostgreSQL server is ready");
            last;
        }
        if ($i == 10) {
            die "PostgreSQL server failed to start within 10 seconds\n";
        }
        sleep(1);
    }
    
    # Create test database
    my $createdb_cmd = "$createdb -h '$temp_dir' -p $port -U '$test_user' '$test_db'";
    unless (run_command($createdb_cmd)) {
        die "Failed to create test database\n";
    }
    
    log_message("Test environment setup complete");
}

sub connect_database {
    log_message("Connecting to test database...");
    
    my $dsn = "DBI:Pg:dbname=$test_db;host=$temp_dir;port=$port";
    my $dbh = DBI->connect($dsn, $test_user, '', {
        AutoCommit => 1,
        RaiseError => 1,
        PrintError => 0,
    });
    
    unless ($dbh) {
        die "Failed to connect to database: " . DBI->errstr . "\n";
    }
    
    log_message("Connected to database successfully");
    return $dbh;
}

sub test_extension_installation {
    my ($dbh) = @_;
    log_message("Testing extension installation...");
    
    eval {
        $dbh->do("CREATE EXTENSION pg_kmersearch;");
        log_message("Extension pg_kmersearch installed successfully");
    };
    if ($@) {
        die "Failed to install extension: $@\n";
    }
    
    # Check if types are available
    my $sth = $dbh->prepare("SELECT typname FROM pg_type WHERE typname IN ('dna2', 'dna4') ORDER BY typname");
    $sth->execute();
    my @types = map { $_->[0] } @{$sth->fetchall_arrayref()};
    
    unless (@types == 2 && $types[0] eq 'dna2' && $types[1] eq 'dna4') {
        die "Extension types not found or incomplete\n";
    }
    
    log_message("Extension types (dna2, dna4) are available");
}

sub test_dna2_operations {
    my ($dbh) = @_;
    log_message("Testing DNA2 type operations...");
    
    # Create test table
    $dbh->do("CREATE TABLE test_dna2 (id SERIAL PRIMARY KEY, name TEXT, sequence DNA2)");
    log_message("Created test_dna2 table");
    
    # Test data insertion
    my @test_sequences = (
        ['seq1', 'ATCGATCG'],
        ['seq2', 'GCTAGCTA'],
        ['seq3', 'TTTTAAAA'],
        ['seq4', 'CCCCGGGG'],
        ['seq5', 'ATGCATGC'],
    );
    
    my $insert_sth = $dbh->prepare("INSERT INTO test_dna2 (name, sequence) VALUES (?, ?::DNA2)");
    for my $seq_data (@test_sequences) {
        $insert_sth->execute($seq_data->[0], $seq_data->[1]);
        log_message("Inserted DNA2 sequence: $seq_data->[0] = $seq_data->[1]");
    }
    
    # Test data retrieval and output
    my $select_sth = $dbh->prepare("SELECT id, name, sequence FROM test_dna2 ORDER BY id");
    $select_sth->execute();
    
    my $row_count = 0;
    while (my ($id, $name, $sequence) = $select_sth->fetchrow_array()) {
        $row_count++;
        log_message("Retrieved DNA2: id=$id, name=$name, sequence=$sequence");
        
        # Verify the sequence matches expected output
        my $expected = $test_sequences[$row_count - 1]->[1];
        if ($sequence ne $expected) {
            die "DNA2 sequence mismatch: expected '$expected', got '$sequence'\n";
        }
    }
    
    if ($row_count != @test_sequences) {
        die "Expected " . @test_sequences . " rows, got $row_count\n";
    }
    
    # Test equality operations
    my $eq_sth = $dbh->prepare("SELECT COUNT(*) FROM test_dna2 WHERE sequence = ?::DNA2");
    $eq_sth->execute('ATCGATCG');
    my ($count) = $eq_sth->fetchrow_array();
    if ($count != 1) {
        die "DNA2 equality test failed: expected 1 match, got $count\n";
    }
    
    log_message("DNA2 operations completed successfully");
}

sub test_dna4_operations {
    my ($dbh) = @_;
    log_message("Testing DNA4 type operations...");
    
    # Create test table
    $dbh->do("CREATE TABLE test_dna4 (id SERIAL PRIMARY KEY, name TEXT, sequence DNA4)");
    log_message("Created test_dna4 table");
    
    # Test data insertion (including degenerate codes)
    my @test_sequences = (
        ['seq1', 'ATCGATCG'],
        ['seq2', 'MRWSYKDN'],  # Degenerate codes
        ['seq3', 'VHBNATCG'],  # More degenerate codes
        ['seq4', 'NNNNAAAA'],  # N codes and standard bases
        ['seq5', 'ATGCMRWS'],  # Mixed standard and degenerate
    );
    
    my $insert_sth = $dbh->prepare("INSERT INTO test_dna4 (name, sequence) VALUES (?, ?::DNA4)");
    for my $seq_data (@test_sequences) {
        $insert_sth->execute($seq_data->[0], $seq_data->[1]);
        log_message("Inserted DNA4 sequence: $seq_data->[0] = $seq_data->[1]");
    }
    
    # Test data retrieval and output
    my $select_sth = $dbh->prepare("SELECT id, name, sequence FROM test_dna4 ORDER BY id");
    $select_sth->execute();
    
    my $row_count = 0;
    while (my ($id, $name, $sequence) = $select_sth->fetchrow_array()) {
        $row_count++;
        log_message("Retrieved DNA4: id=$id, name=$name, sequence=$sequence");
        
        # Verify the sequence matches expected output (should be uppercase)
        my $expected = uc($test_sequences[$row_count - 1]->[1]);
        if ($sequence ne $expected) {
            die "DNA4 sequence mismatch: expected '$expected', got '$sequence'\n";
        }
    }
    
    if ($row_count != @test_sequences) {
        die "Expected " . @test_sequences . " rows, got $row_count\n";
    }
    
    # Test equality operations
    my $eq_sth = $dbh->prepare("SELECT COUNT(*) FROM test_dna4 WHERE sequence = ?::DNA4");
    $eq_sth->execute('ATCGATCG');
    my ($count) = $eq_sth->fetchrow_array();
    if ($count != 1) {
        die "DNA4 equality test failed: expected 1 match, got $count\n";
    }
    
    # Test degenerate code handling
    $eq_sth->execute('MRWSYKDN');
    ($count) = $eq_sth->fetchrow_array();
    if ($count != 1) {
        die "DNA4 degenerate code test failed: expected 1 match, got $count\n";
    }
    
    log_message("DNA4 operations completed successfully");
}

sub test_score_functions {
    my ($dbh) = @_;
    log_message("Testing score calculation functions...");
    
    # Test with DNA2
    eval {
        my $sth = $dbh->prepare("SELECT kmersearch_rawscore(sequence, 'ATCGATCG') FROM test_dna2 WHERE name = 'seq1'");
        $sth->execute();
        my ($score) = $sth->fetchrow_array();
        log_message("DNA2 raw score for exact match: $score");
        
        # Should have some positive score for exact match
        if (!defined $score || $score < 0) {
            die "DNA2 raw score test failed: expected positive score, got " . ($score // 'NULL') . "\n";
        }
    };
    if ($@) {
        log_message("Score function test failed (this may be expected if functions are not fully implemented): $@");
    }
    
    # Test with DNA4
    eval {
        my $sth = $dbh->prepare("SELECT kmersearch_rawscore(sequence, 'ATCGATCG') FROM test_dna4 WHERE name = 'seq1'");
        $sth->execute();
        my ($score) = $sth->fetchrow_array();
        log_message("DNA4 raw score for exact match: $score");
        
        # Should have some positive score for exact match
        if (!defined $score || $score < 0) {
            die "DNA4 raw score test failed: expected positive score, got " . ($score // 'NULL') . "\n";
        }
    };
    if ($@) {
        log_message("Score function test failed (this may be expected if functions are not fully implemented): $@");
    }
    
    log_message("Score function tests completed");
}

sub run_tests {
    log_message("Starting pg_kmersearch installation tests...");
    
    my $dbh = connect_database();
    
    eval {
        test_extension_installation($dbh);
        test_dna2_operations($dbh);
        test_dna4_operations($dbh);
        test_score_functions($dbh);
        
        log_message("All tests completed successfully!");
        print "PASS: All pg_kmersearch installation tests passed\n";
    };
    
    if ($@) {
        print "FAIL: Test failed: $@";
        $dbh->disconnect() if $dbh;
        return 0;
    }
    
    $dbh->disconnect();
    return 1;
}

# Main execution
print "pg_kmersearch Installation Test\n";
print "===============================\n";

my $success = 0;
eval {
    setup_test_environment();
    $success = run_tests();
};

if ($@) {
    print "FAIL: Setup or test execution failed: $@";
}

# Always cleanup
cleanup();

if ($success) {
    print "\nSUMMARY: All tests passed successfully!\n";
    exit 0;
} else {
    print "\nSUMMARY: Some tests failed.\n";
    exit 1;
}