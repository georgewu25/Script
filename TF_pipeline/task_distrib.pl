#!/bin/perl
# use this script to distribute evenly among tasks
# to prepare for array job

use strict;
use POSIX qw(ceil);

my $USAGE=<<EOF;
Usage: perl task_distrib.pl task_list total_task task_id 
--task_list: the file contains the list of items to be grouped
--total_task: total task requested for job
--task_id: current task id
EOF

if($#ARGV<2) {
    print "$USAGE\n";
}
my ($task_list, $total_task, $task_id)=@ARGV;
my @task_list_arr=();
# get the total fastq file list:
if( -d $task_list ) { # task_list is given using fasta folder
    opendir (DIR, $task_list) or die "Failed to open directory\n";
    @task_list_arr = grep { /fastq/ } readdir DIR;
    closedir (DIR);
}
elsif (-e $task_list) { # task list is given in a file
    open IN, "<$task_list" || die "Error in opening $task_list!";
    # while (my $line=<IN>) {
    # 	chomp $line;
    # 	my ($sid, $sample_name)=split(",", $line);
    # 	push @task_list_arr, $sample_name;
    # }
    @task_list_arr=<IN>;
    close IN;
}

my $total_size=$#task_list_arr + 1;
my $batch_size = ceil($total_size/$total_task);
# fix the bug:
#my $last_batch = $total_size% $total_task;
my $last_batch = ($total_size-1) % $total_task+1;
my $adj_batch=$batch_size-$last_batch;
my @batch=();

if($task_id<=$last_batch) { # full size
  @batch = splice(@task_list_arr, $batch_size*($task_id-1), $batch_size);
}
else {
  @batch = splice(@task_list_arr, ($batch_size-1)*($task_id-1)+$last_batch, $batch_size-1);  
}

print join "\n", @batch;


