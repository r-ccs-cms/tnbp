#!/usr/bin/env perl
# Usage: ./getdata.p [input_file]
use strict;
use warnings;

my $infile = $ARGV[0] // "output.dat";

open my $fh, "<", $infile or die "Cannot open $infile: $!";

my (@x, @z, @zz);

while (my $line = <$fh>) {

    # X成分
    if ($line =~ /measurement of X\b/) {
        my ($ts)   = $line =~ /time step\s+(\d+)/;
        my ($site) = $line =~ /site\s+(\d+)/;
        push @x, [$ts // -1, $site // -1, $line] if defined $ts && defined $site;
    }

    # Z成分 (スペース入りを区別)
    elsif ($line =~ /measurement of Z\s/) {
        my ($ts)   = $line =~ /time step\s+(\d+)/;
        my ($site) = $line =~ /site\s+(\d+)/;
        push @z, [$ts // -1, $site // -1, $line] if defined $ts && defined $site;
    }

    # ZZ相関
    elsif ($line =~ /measurement of ZZ/) {
        my ($ts) = $line =~ /time step\s+(\d+)/;
        my ($i,$j) = $line =~ /sites\s*\((\d+),\s*(\d+)\)/;
        push @zz, [$ts // -1, $i // -1, $line, ($j // -1)] if defined $ts && defined $i;
    }
}
close $fh;

# --- ソート ---
# X: time step -> site
@x = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @x;

# Z: time step -> site
@z = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @z;

# ZZ: time step -> i -> j
@zz = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] || $a->[3] <=> $b->[3] } @zz;

# --- 出力 ---
open my $ox, ">", "data_meas_x.txt" or die $!;
print $ox $_->[2] for @x;
close $ox;

open my $oz, ">", "data_meas_z.txt" or die $!;
print $oz $_->[2] for @z;
close $oz;

open my $oj, ">", "data_meas_j.txt" or die $!;
print $oj $_->[2] for @zz;
close $oj;

print "Wrote: data_meas_x.txt, data_meas_z.txt, data_meas_j.txt\n";
