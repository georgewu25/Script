use strict;

my $fc_clean="featureCounts_clean.txt";
my $fc_reorder="featureCounts_reorder.txt";
my $sample_order="../FASTQ/sample_list.csv";
my @final_cols;

open SO, "<$sample_order";
while (<SO>) {
    chomp $_;
    push @final_cols, $_ if $_ ne ""; 
}
close SO;
open FCR, ">$fc_reorder";
open FC, "<$fc_clean";
my $line = <FC>;
my %fc = ();
chomp $line;
my @colnames = split("\t", $line);

print FCR join "\t", $colnames[0], @final_cols;
print FCR "\n";

while ($line=<FC>) {
    my @fc_val = split("\t", $line);
    @fc{@colnames}=@fc_val;
    print FCR join "\t", $fc{Geneid}, @fc{@final_cols};
    print FCR "\n";
}
close FC;
close FCR;

print "Done, the result is in $fc_reorder\n";
