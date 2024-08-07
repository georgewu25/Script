#!/bin/perl

use strict;

# this is to clean up the featurecount result file
# 1) remove column 2-6
# 2) remove first line
# 3) clean up the headline with sample name only as column name

my $USAGE= <<EOF;

Usage: perl clean_fc_headline.pl featurecount_output_file
command line argument: 
'featurecount_output_file' - the featurecount output file from 'featurecount' program
EOF

if($#ARGV<1) {
    print "$USAGE\n";
    exit;
}
my ($fq_orig_dir, $fq_new_dir)=@ARGV; 
