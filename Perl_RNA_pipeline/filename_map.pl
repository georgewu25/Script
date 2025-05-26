#!/bin/perl

use strict;

# this is to rewrite the sample fastq map
# by simply mapping to the desired order 
# with the provided sample order file
# the sample order file will be named as:
# {PROJ}_ordered_sample.txt
# originally, this file will have the format: 
# sample1_basename
# sample2_basename

my $USAGE= <<EOF;

Usage: perl filename_map_apoe22.pl fq_orig_dir fq_new_dir


EOF

if($#ARGV<1) {
    print "$USAGE\n";
    exit;
}

my ($fq_orig_dir, $fq_new_dir)=@ARGV; 

# get the total fastq file list in the $fq_orig_dir:


opendir(my $dh, $fq_orig_dir) || die "Can't opendir $fq_orig_dir: $!";
my @fq_orig_list = grep { /(R1|R2)_001\.fastq\.gz$/ && -f "$fq_orig_dir/$_" } readdir($dh);
closedir $dh;

system("mkdir -p $fq_new_dir");
open SMAP, ">$fq_new_dir/sample_map.csv";
print SMAP join ",", "fq_id", "sample_basename", 'read', "lane", "orig_fname", "slink_fname";
print SMAP "\n";
open SERR, ">$fq_new_dir/sample_err.csv";
print SERR join ",", "fq_id", "sample_basename", 'read', "lane", "orig_fname", "slink_conflict", "exit_code", "errmsg";
print SERR "\n";
my $id=0; # id for each fastq file
my %sample_list=(); # put together all the sample names for downstream analysis
my %sample_name_list=(); # put together each unique sequence file name, remove region info for alignment
my $bid=0; # sid for each unique sample base
my $sid=0; # sid for each unique sample name, embeded with treatment info
my $sample_basename;
my $sample_name;
# my $r1_count=0;
# my $r2_count=0;
my $ftype='fastq';
my $gz='gz';

foreach my $fq (sort @fq_orig_list) {
    my ($ID, $read, $lane, $filetype, $gz);  # Declare variables
    
    my @items = split("[_.]", $fq);

    $ID = $items[0];
    $read = $items[-4];
    $lane = $items[-3];
    $filetype = $items[-2];
    $gz = $items[-1];

    $sample_basename = $ID;

    $sample_list{$sample_basename}{$read}{$lane}{ORIG_FILENAME} = $fq; # Defines ORIG_FILENAME

    print join "::", $ID, $read, $lane, $filetype, $gz;
    print "\n";
}



foreach my $sample_basename (sort keys %sample_list){
  $bid++; 
  foreach my $read (sort keys %{$sample_list{$sample_basename}}) {
    foreach my $lane (sort keys %{$sample_list{$sample_basename}{$read}}) {
    $id++; # give each fastq file an id
	  $sample_list{$sample_basename}{$read}{$lane}{BASE_ID}=$bid; # sample root(base) id
	  $sample_list{$sample_basename}{$read}{$lane}{NEW_FILENAME}=(join "_", $id, $sample_basename, $read, $lane) . ".$ftype" . ".$gz";
 
	  print "NEW FILENAME: " . $sample_list{$sample_basename}{$read}{$lane}{NEW_FILENAME};
          print "\n";
          $sample_list{$sample_basename}{$read}{$lane}{FQ_ID}=$id;

    my $rtn=system("ln -s $fq_orig_dir/$sample_list{$sample_basename}{$read}{$lane}{ORIG_FILENAME} $fq_new_dir/$sample_list{$sample_basename}{$read}{$lane}{NEW_FILENAME}");
    
	    if($rtn == 0) {
	    print SMAP join ",", 
	    $sample_list{$sample_basename}{$read}{$lane}{FQ_ID},
        $sample_basename,
	    $sample_list{$sample_basename}{$read}{$lane}{BASE_ID},
	    $lane,
	    $sample_list{$sample_basename}{$read}{$lane}{ORIG_FILENAME},
  	    $sample_list{$sample_basename}{$read}{$lane}{NEW_FILENAME};
            print SMAP "\n";
  	  } 
          else { # maybe already linked
	    print SERR join ",", 
	    $sample_list{$sample_basename}{$read}{$lane}{FQ_ID},
	    $sample_basename,
	    $sample_list{$sample_basename}{$read}{$lane}{BASE_ID},
	    $lane,
	    $sample_list{$sample_basename}{$read}{$lane}{ORIG_FILENAME},
  	    $sample_list{$sample_basename}{$read}{$lane}{NEW_FILENAME},
 	    $?>>8, $!;
	    print SERR "\n"; 	
          } # end if-else
        } # end lane
        } #end read
}  # basename

close SMAP;
close SERR;

open STAT, ">$fq_new_dir/sample_file_statisitics.txt";
print STAT "total fastq files processed: $id";
print STAT "\n";
print STAT "total sample count: $bid";
print STAT "\n";
close STAT;

# output the unique sample base names:
open SBLIST, ">$fq_new_dir/samplebase_list.csv";
for my $s (sort keys %sample_list) {
    print "sample base: $s\n";
    print SBLIST $s;
    print SBLIST "\n";
}
close SBLIST;

# output the unique sequence names: 
open SLIST, ">$fq_new_dir/sample_treatment_list.csv";
for my $st (sort {$a <=> $b } keys %sample_list) {
    print SLIST $st;
    print SLIST "\n";
}
close SLIST;

