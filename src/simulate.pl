#! /usr/bin/env perl

use strict;
use Math::Random;

my $seq_len    = 2000;
my $max_iter   = 20000;
my $insert_len = shift @ARGV;
my $read_len   = 100;

my @counts;

foreach ( 1 .. $max_iter ) {

    # Generate a valid pair of reads

    my $left;
    my $right;

    while (1) {
        $left = int( rand($seq_len) );
        next if $left + $read_len >= $seq_len;

        $right = $left + Math::Random::random_poisson( 1, $insert_len );
        next if $right + $read_len >= $seq_len;

        last;
    }

    foreach my $pos ( $left .. $right + $read_len ) {
        $counts[$pos] += 1;
    }
}

foreach my $pos ( 0 .. $#counts ) {
    print join( "\t", $pos, $counts[$pos] ) . "\n";
}
