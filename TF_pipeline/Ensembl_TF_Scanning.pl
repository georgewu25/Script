#!/usr/bin/perl

use strict;
use warnings;
use File::Path qw(make_path);

sub Ensembl_setup {
    my ($query_name, $brain_region, $cell_type, $ADNC, $input_gene_txt, $input_gene_fasta, $background_gene_txt, $background_gene_fasta) = @_;
    
    # Define main directory for output
    my $output_main_dir = "/Project/SEA-AD_ATAC/$brain_region/$cell_type/$ADNC/";

    #### MORA setup
    # MORA requires the fasta files inside its local dir, which it will generate con files based on the fa files
    # For file organization, we save the input files in the input folder, copy them to the local dir, and remove once finished.
    my $mora_main_dir = "/Project/TF_Analysis/MORA/";

    # Read template qsub
    my $mora_template_qsub_path = "$mora_main_dir/MORA_qsub/Template.qsub";
    my @mora_template_qsub;
    {
        open my $fh, '<', $mora_template_qsub_path or die "Could not open '$mora_template_qsub_path': $!\n";
        @mora_template_qsub = <$fh>;
        close $fh;
    }
    
    ## MORA requires the input files to be present in the same dir as its scripts, so I make copies and delete afterwards
    my $mora_script_dir

    # Replace placeholders
    @mora_template_qsub = map { s/job_placeholder/$query_name/gr } @mora_template_qsub;
    @mora_template_qsub = map { s/input_placeholder/"$input_gene_fasta"/gr } @mora_template_qsub;
    @mora_template_qsub = map { s/background_placeholder/"$background_gene_fasta"/gr } @mora_template_qsub;
    @mora_template_qsub = map { s/output_dir_placeholder/"$output_main_dir\/MORA_output/"/gr } @mora_template_qsub;
    

    my $mora_final_qsub = "$mora_main_dir/MORA_qsub/$query_name.qsub";
    {
      open my $fh, '>', $mora_final_qsub or die "Could not open '$mora_final_qsub': $!\n";
      print $fh join("\n", @mora_template_qsub);
      close $fh;
    }

    chmod 0755, $mora_final_qsub or die "Failed to make MORA final script executable: $!\n";
    system("qsub $mora_final_qsub") == 0 or die "MORA qsub failed: $!\n";
    

    #### AME setup
    # AME also requires gene sequence fasta files, which can be the same as those for MORA
    my $ame_main_dir = "/Project/TF_Analysis/AME/";

    # Read template qsub
    my $ame_template_qsub_path = "$ame_main_dir/Template.qsub";
    my @ame_template_qsub;
    {
        open my $fh, '<', $ame_template_qsub_path or die "Could not open '$ame_template_qsub_path': $!\n";
        @ame_template_qsub = <$fh>;
        close $fh;
    }

    # Replace placeholders
    @ame_template_qsub = map { s/job_placeholder/$query_name/gr } @ame_template_qsub;
    @ame_template_qsub = map { s/input_placeholder/"$input_gene_fasta"/gr } @ame_template_qsub;
    @ame_template_qsub = map { s/output_dir_placeholder/"$output_main_dir\/AME_output/"/gr } @ame_template_qsub;

    my $ame_final_qsub = "$ame_main_dir/$query_name.qsub";
    {
      open my $fh, '>', $ame_final_qsub or die "Could not open '$ame_final_qsub': $!\n";
      print $fh join("\n", @ame_template_qsub);
      close $fh;
    }
    
    chmod 0755, $ame_final_qsub or die "Failed to make AME final script executable: $!\n";
    system("qsub $ame_final_qsub") == 0 or die "AME qsub failed: $!\n";
    
    
    #### HOMER Setup (HOMER takes gene names as inputs)
    my $homer_main_dir = "/Project/TF_Analysis/HOMER/";
    
     # Read template qsub
    my $homer_template_path = "$homer_main_dir/Template.sh";
    my @homer_template;
    {
        open my $fh, '<', $homer_template_path or die "Could not open '$homer_template_path': $!\n";
        @homer_template = <$fh>;
        close $fh;
    }

    # Replace placeholders
    @homer_template = map { s/input_placeholder/"$input_gene_txt"/gr } @homer_template;
    @homer_template = map { s/output_dir_placeholder/"$output_main_dir\/HOMER_output/"/gr } @homer_template;
    
    my $homer_final = "$homer_main_dir/$query_name.sh";
    {
      open my $fh, '>', $homer_final or die "Could not open '$homer_final': $!\n";
      print $fh join("\n", @homer_template);
      close $fh;
    }
    
    chmod 0755, $homer_final or die "Failed to make HOMER final script executable: $!\n";
    system("bash -c \"$homer_final\"") == 0 or die "HOMER Submission Fail: $!\n";
    


    #### Bart2 Setup (Bart2 takes gene names as inputs like HOMER)
    my $bart_main_dir = "/Project/TF_Analysis/Bart2/";
    
    # Read template qsub
    my $bart_template_path = "$bart_main_dir/Template.sh";
    my @bart_template;
    {
        open my $fh, '<', $bart_template_path or die "Could not open '$bart_template_path': $!\n";
        @bart_template = <$fh>;
        close $fh;
    }
    
    # Replace placeholders
    @bart_template = map { s/input_placeholder/"$input_gene_txt"/gr } @bart_template;
    @bart_template = map { s/output_dir_placeholder/"$output_main_dir\/Bart2_output/"/gr } @bart_template;
    
    my $bart_final = "$bart_main_dir/$query_name.sh";
    {
      open my $fh, '>', $bart_final or die "Could not open '$bart_final': $!\n";
      print $fh join("\n", @bart_template);
      close $fh;
    }
    
    chmod 0755, $bart_final or die "Failed to make Bart2 final script executable: $!\n";
    system("bash -c \"$bart_final\"") == 0 or die "Bart2 Submission Fail: $!\n";
    
    
    ####Lisa2 setup (Lisa takes gene names from both query and background as inputs)
    my $lisa_main_dir = "/Project/TF_Analysis/Lisa/";
    
    # Read template qsub
    my $lisa_template_path = "$lisa_main_dir/Template.sh";
    my @lisa_template;
    {
        open my $fh, '<', $lisa_template_path or die "Could not open '$lisa_template_path': $!\n";
        @lisa_template = <$fh>;
        close $fh;
    }
    
    # Replace placeholders
    @lisa_template = map { s/input_placeholder/"$input_gene_txt"/gr } @lisa_template;
    @lisa_template = map { s/background_placeholder/"$background_gene_txt"/gr } @lisa_template;
    @lisa_template = map { s/output_prefix_placeholder/"$output_main_dir\/Lisa_output"/gr } @lisa_template;
    
    my $lisa_final = "$lisa_main_dir/$query_name.sh";
    {
      open my $fh, '>', $lisa_final or die "Could not open '$lisa_final': $!\n";
      print $fh join("\n", @lisa_template);
      close $fh;
    }
    
    chmod 0755, $lisa_final or die "Failed to make Lisa2 final script executable: $!\n";
    system("bash -c \"$lisa_final\"") == 0 or die "Lisa2 Submission Fail: $!\n";
    
    #### Pscan Setup (Pscan has a web-browser and I have not downloaded the source code, it also takes gene Refseq IDs as inputs)
}