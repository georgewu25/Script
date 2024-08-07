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
# ...
# 1. Study, (Astro or Micro/A or M)
# 2. Individual, (TCW1 or TCW2)
# 3. Genotype (33(A/B/C) and 44(a/b/c))
# 4. Region (R1/R2)
# 5. Lane (001)
# 4. Treatment Condition (My/control)
# 5. Treatment hours (4h, 20h)
# sample name: 
# A-TCW2-33C-20h_R2_001.fastq.gz
# 
my $USAGE= <<EOF;

Usage: perl filename_map_all_3.pl fq_orig_dir fq_new_dir ordered_sample

It will
1) create the softlink for each fastq.gz file under 
fq_orig_dir, using the name as {sample_id}_{norm_sample}_READ{1|2}.fastq.gz
2) generate a sample_map.txt file contains all sample id 
maps;
3) report any erroreous fastq.gz that can\'t be softlinked and 
put in sample_err.txt
4) generate a text file \'sample_list.txt\' to list all unique
sample id in the directory. 
EOF

if($#ARGV<2) {
    print "$USAGE\n";
    exit;
}
# test with /projectnb/tcwlab/RawRNAseq/APOE22/00_fastq_100
my ($fq_orig_dir, $fq_new_dir, $sample_file)=@ARGV; 
# get the total fastq file list in the $fq_orig_dir:
open(my $dh, $sample_file) || die "Can't opendir $sample_file: $!";
close $dh;
my $sample_str = do {
    local $/ = undef;
    open my $fh, "<", $sample_file
        or die "could not open $sample_file: $!";
    <$fh>;
};
my %samples=split(",|\n", $sample_str);

opendir(my $dh, $fq_orig_dir) || die "Can't opendir $fq_orig_dir: $!";
my @fq_orig_list = grep { /\.fastq\.gz/ && -f "$fq_orig_dir/$_" } readdir($dh);
closedir $dh;

system("mkdir -p $fq_new_dir");
open SMAP, ">$fq_new_dir/sample_map.csv";
print SMAP join ",", "fq_id", "sample_basename", "base_id", "sample_name", "sample_id", "study", "indv", "genotype", "treatmc", "treatmh", "treatmch", "region", "lane", "orig_fname", "slink_fname";
print SMAP "\n";
open SERR, ">$fq_new_dir/sample_err.csv";
print SERR join ",", "fq_id", "sample_basename", "base_id", "sample_name", "sample_id", "study", "indv", "genotype", "treatmc", "treatmh", "treatmch", "region", "lane", "orig_fname", "slink_conflict", "exit_code", "errmsg";
print SERR "\n";
my $id=0; # id for each fastq file
my %sample_list=(); # put together all the sample names for downstream analysis
my %sample_name_list=(); # put together each unique sequence file name, remove region info for alignment
my $bid=0; # sid for each unique sample base
my $sid=0; # sid for each unique sample name, embeded with treatment info
my $sample_basename;
my $sample_name;
#my $tid=0; 
# my $r1_count=0;
# my $r2_count=0;
my $ftype='fastq';
my $gz='gz';

foreach my $fq (sort @fq_orig_list) {
    $sample_basename="";
  # APOE_Mye fastq name
#  file name example: M-TCW1-44a-My20h_R1_001.fastq.gz #ATGTCA_M_TCW1_44a-READ1.fastq.gz
#  my ($study, $indv, $genotype, $treatmch, $region, $lane, $file_type, $gz) = split("[-_.]", $fq);
    # APOE22 fastq name: M2-22-1_R1_001.fastq.gz (study, genotype, indv, region, lane
    #                    or M2-22-4h-1_R1_001.fastq.gz  (study,genotype, tmch, indv. 
    #  my ($study, $genotype, $treatmch, $indv, $region, $lane, $file_type, $gz) = split("[-_.]", $fq);
    my ($study, $genotype, $treatmch, $indv, $region, $lane, $file_type, $gz);
    my ($treatmc, $treatmh);
    my @items = split("[-_.]", $fq);
    if ($#items==7) {
	$study=$items[0];
	$genotype=$items[1];
	$treatmch=$items[2];
	$indv=$items[3];
	$region=$items[4];
	$lane=$items[5];
	$file_type=$items[6];
	$gz=$items[7];
	($treatmc, $treatmh) = $treatmch=~/([My]*)(\d+)h/;
	$treatmc="Ctrl" if $treatmc eq "";
    }
    elsif($#items==6) {
	$study=$items[0];
	$genotype=$items[1];
	$treatmch="";
	$treatmc="Ctrl";
	$treatmh="";
	$indv=$items[2];
	$region=$items[3];
	$lane=$items[4];
	$file_type=$items[5];
	$gz=$items[6];	
    }
    else {
	print SERR "$fq is not part of expected file name format. Please deal with it off line. \n"
    }
  my $studyindv = join "_", $study, $indv;  
  $sample_basename=join "_", $study, $genotype;
    if($treatmch ne "") {
	$sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{SAMPLE_NAME}=join "_", $sample_basename, $treatmch, $indv;
    }
    else {
	$sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{SAMPLE_NAME}=join "_", $sample_basename, $indv;	
    }
  $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{TREATMCH}=$treatmch;
  $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{GENOTYPE}=$genotype;
  $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{STUDY}=$study;
  $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{STUDYINDV}=$studyindv;
  $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{ORIG_FILENAME}=$fq;

    print join "::", $study, $genotype, $treatmch, $indv, $region, $lane, $file_type, $gz;
    print "\n";
   
}

foreach my $sample_basename (sort keys %sample_list){
  $bid++; 
  foreach my $treatmc (sort keys %{$sample_list{$sample_basename}}) {
    foreach my $treatmh (sort {$a<=>$b} keys  %{$sample_list{$sample_basename}{$treatmc}}) {
      foreach my $indv (sort keys %{$sample_list{$sample_basename}{$treatmc}{$treatmh}}) {
	$sid++;
        foreach my $region (sort keys %{$sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}}) {
          foreach my $lane (sort keys %{$sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}}) {
	  $id++; # give each fastq file an id
	  $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{BASE_ID}=$bid; # sample root(base) id
	  $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{SAMPLE_ID}=$sid; # same sample can have multiple fq's
	  $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{NEW_FILENAME}
                   =(join "_", $id, $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{SAMPLE_NAME}, $region, $lane) . ".$ftype" . ".$gz"; 
	  print "NEW FILENAME: " . $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{NEW_FILENAME};
          print "\n";
          $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{FQ_ID}=$id;
	  # record sample name in right order:
          $sample_name_list{$sid} = $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{SAMPLE_NAME};

	  my $rtn=system("ln -s $fq_orig_dir/$sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{ORIG_FILENAME} $fq_new_dir/$sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{NEW_FILENAME}");
	  if($rtn == 0) {
	    print SMAP join ",", 
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{FQ_ID},
            $sample_basename,
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{BASE_ID},
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{SAMPLE_NAME},
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{SAMPLE_ID},
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{STUDY},
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{INDIVIDUAL},
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{GENOTYPE},
	    $treatmc,
	    $treatmh,
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{TREATMCH},
            $region,
	    $lane,
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{ORIG_FILENAME},
  	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{NEW_FILENAME};
            print SMAP "\n";
  	  } 
          else { # maybe already linked
	    print SERR join ",", 
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{FQ_ID},
	    $sample_basename,
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{BASE_ID},
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{SAMPLE_NAME},
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{SAMPLE_ID},
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{STUDY},
	    $indv,
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{GENOTYPE},
	    $treatmc,
	    $treatmh,
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{TREATMCH},
            $region,
	    $lane,
	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{ORIG_FILENAME},
  	    $sample_list{$sample_basename}{$treatmc}{$treatmh}{$indv}{$region}{$lane}{NEW_FILENAME},
 	    $?>>8, $!;
	    print SERR "\n"; 	
          } # end if-else
        } # end lane
      } # end region
      } # end indv
    } # end mh	
  } # end mc
}  # basename

close SMAP;
close SERR;

open STAT, ">$fq_new_dir/sample_file_statisitics.txt";
print STAT "total fastq files processed: $id";
print STAT "\n";
print STAT "total sample count: $sid";
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
open SLIST, ">$fq_new_dir/sample_name_list.csv";
for my $st (sort {$a <=> $b } keys %sample_name_list) {
    print SLIST join ",", $st, $sample_name_list{$st};
    print SLIST "\n";
}
close SLIST;
