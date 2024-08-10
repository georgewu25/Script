#!/usr/bin/env perl 
use strict;
my $USAGE= <<EOF;

Usage: perl generate_proj_script_w_config.pl proj_config_file
command line argument: 
'proj_config_file' - the configuration file contains all the meta data for the pipeline 


if($#ARGV<0) {
    print "$USAGE\n";
    exit;
}

my %proj_list = ();
my $curproj="";

# read input from config file: 
my $config=shift @ARGV;
open CFG, "<$config";
while(my $line=<CFG>) {
    next if $line=~/^#/; # skip comment line
    if($line=~/\[(.+)\]/) {
      $curproj = $1;
      next;
    }
    if($line=~/(.+):(.+)/) {
    	$proj_list{$curproj}{$1}=$2;
	    next;
    }   
}
close CFG;



# get parameters: 

my $QLNK_TEMPLATE="RNAseq_fastq_link.qsub";
my $QC_TEMPLATE="RNAseq_fastqc.qsub";
my $MQC_TEMPLATE="RNAseq_multiqc.qsub";
my $ALN_TEMPLATE="RNAseq_aln.qsub";
my $MQC_ALN_TEMPLATE="RNAseq_multiqc_aln.qsub";
my $IDX_BAM_TEMPLATE="RNAseq_idx_bam.qsub";
my $FC_TEMPLATE="RNAseq_fc.qsub";
my $DS2_TEMPLATE="RNAseq_deseq2.qsub";
my $DS2_R_TEMPLATE="RNAseq_deseq2.R";
my $GSEA_TEMPLATE="RNAseq_gsea.qsub";
my $GSEA_R_TEMPLATE="RNAseq_gsea.R";
my $PIPE_TEMPLATE="RNAseq_pipeline.sh";

my $script_dir;
my $QLNK_QSUB;
my $QC_QSUB;
my $MQC_QSUB;
my $ALN_QSUB;
my $MQC_ALN_QSUB;
my $IDX_BAM_QSUB;
my $FC_QSUB;
my $DS2_QSUB;
my $DS2_R;
my $GSEA_QSUB;
my $GSEA_R;
my $PIPE_SCRIPT;

=pod


foreach my $proj (keys %proj_list) {
    $script_dir = $proj_list{$proj}{PROJ_DIR} . "/scripts/";
    system("mkdir -p $script_dir");
    $QLNK_QSUB= $script_dir . "RNAseq_fastq_link_" . $proj . ".qsub";
    $QC_QSUB= $script_dir . "RNAseq_fastqc_" . $proj . ".qsub";
    $MQC_QSUB= $script_dir . "RNAseq_multiqc_" . $proj . ".qsub";
    $ALN_QSUB= $script_dir . "RNAseq_aln_" . $proj . ".qsub";
    $MQC_ALN_QSUB= $script_dir . "RNAseq_multiqc_aln_" . $proj . ".qsub";
    $IDX_BAM_QSUB= $script_dir . "RNAseq_idx_bam_" . $proj . ".qsub";
    $FC_QSUB= $script_dir . "RNAseq_fc_" . $proj . ".qsub";
    $DS2_QSUB= $script_dir . "RNAseq_deseq2_" . $proj . ".qsub";
    $DS2_R= $script_dir . "RNAseq_deseq2_" . $proj . ".R";
    $GSEA_QSUB= $script_dir . "RNAseq_gsea_" . $proj . ".qsub";    
    $GSEA_R= $script_dir . "RNAseq_gsea_" . $proj . ".R";    
    $PIPE_SCRIPT= $script_dir  . "RNAseq_pipeline_". $proj . ".sh";

# step 1: generate softlink bash: 
    open IN, "<$QLNK_TEMPLATE";
    my $content = do { local $/; <IN> };
    close IN;
    $str_to_replace{COMMON_SCRIPT_DIR} = $common_script_dir;
    $str_to_replace{PROJECTNAME} = $proj;
    $str_to_replace{PROJECT_RAWDATA_DIR} = $proj_list{$proj}{PROJ_RAWDIR};
    $str_to_replace{PROJECT_DIR} = $proj_list{$proj}{PROJ_DIR};
    $str_to_replace{FILEMAP_SCRIPT} = $proj_list{$proj}{PROJ_MAP_SCRIPT};
    for my $k (keys %str_to_replace) {
    	$content =~ s/\{$k\}/$str_to_replace{$k}/g;
    }

    open OUT, ">$QLNK_QSUB";
    print OUT $content;
    close OUT;

# prepare all the subdirectories to call qsub: 
    system("mkdir -p $proj_list{$proj}{PROJ_DIR}/qlog/fqlink");

# step 2: generate fastqc qsub: 

    open IN, "<$QC_TEMPLATE";
    my $content = do { local $/; <IN> };
    close IN;
    $str_to_replace{COMMON_SCRIPT_DIR} = $common_script_dir;
    $str_to_replace{PROJECTNAME} = $proj;
    $str_to_replace{PROJECT_RAWDATA_DIR} = $proj_list{$proj}{PROJ_RAWDIR};
    $str_to_replace{PROJECT_DIR} = $proj_list{$proj}{PROJ_DIR};
    $str_to_replace{FILEMAP_SCRIPT} = $proj_list{$proj}{PROJ_MAP_SCRIPT};
    $str_to_replace{SGE_TOTAL_TASK} = $proj_list{$proj}{SGE_TOTAL_TASK};
    for my $k (keys %str_to_replace) {
    	$content =~ s/\{$k\}/$str_to_replace{$k}/g;
    }

    open OUT, ">$QC_QSUB";
    print OUT $content;
    close OUT;

# prepare all the subdirectories to call qsub: 
    system("mkdir -p $proj_list{$proj}{PROJ_DIR}/qlog/fqc");
    #system("qsub $QC_QSUB")


# step 3: multiqc script: 
    %str_to_replace=();

    open IN, "<$MQC_TEMPLATE";
    $content = do { local $/; <IN> };
    close IN;
    $str_to_replace{COMMON_SCRIPT_DIR} = $common_script_dir;
    $str_to_replace{PROJECT_DIR} = $proj_list{$proj}{PROJ_DIR};
    $str_to_replace{PROJECTNAME} = $proj;
    for my $k (keys %str_to_replace) {
    	$content =~ s/\{$k\}/$str_to_replace{$k}/g;
    }

    open OUT, ">$MQC_QSUB";
    print OUT $content;
    close OUT;

    system("mkdir -p $proj_list{$proj}{PROJ_DIR}/qlog/mqc");

# step 4: alignment script: 
    %str_to_replace=();

    open IN, "<$ALN_TEMPLATE";
    $content = do { local $/; <IN> };
    close IN;
    $str_to_replace{COMMON_SCRIPT_DIR} = $common_script_dir;
    $str_to_replace{PROJECT_DIR} = $proj_list{$proj}{PROJ_DIR};
    $str_to_replace{PROJECTNAME} = $proj;
    $str_to_replace{SGE_TOTAL_TASK} = $proj_list{$proj}{SGE_TOTAL_TASK};
    for my $k (keys %str_to_replace) {
      	$content =~ s/\{$k\}/$str_to_replace{$k}/g;
    }
    
    open OUT, ">$ALN_QSUB";
    print OUT $content;
    close OUT;

    system("mkdir -p $proj_list{$proj}{PROJ_DIR}/qlog/aln");

# step 5: multiqc alignment result script: 
    %str_to_replace=();

    open IN, "<$MQC_ALN_TEMPLATE";
    $content = do { local $/; <IN> };
    close IN;
    $str_to_replace{COMMON_SCRIPT_DIR} = $common_script_dir;
    $str_to_replace{PROJECT_DIR} = $proj_list{$proj}{PROJ_DIR};
    $str_to_replace{PROJECTNAME} = $proj;
    for my $k (keys %str_to_replace) {
      	$content =~ s/\{$k\}/$str_to_replace{$k}/g;
    }

    open OUT, ">$MQC_ALN_QSUB";
    print OUT $content;
    close OUT;

    system("mkdir -p $proj_list{$proj}{PROJ_DIR}/qlog/mqc_aln");

# step 6, indexing bam result: 
    %str_to_replace=();

    open IN, "<$IDX_BAM_TEMPLATE";
    $content = do { local $/; <IN> };
    close IN;
    $str_to_replace{COMMON_SCRIPT_DIR} = $common_script_dir;
    $str_to_replace{PROJECT_DIR} = $proj_list{$proj}{PROJ_DIR};
    $str_to_replace{PROJECTNAME} = $proj;
    $str_to_replace{SGE_TOTAL_TASK} = $proj_list{$proj}{SGE_TOTAL_TASK};
    for my $k (keys %str_to_replace) {
      	$content =~ s/\{$k\}/$str_to_replace{$k}/g;
    }

    open OUT, ">$IDX_BAM_QSUB";
    print OUT $content;
    close OUT;

    system("mkdir -p $proj_list{$proj}{PROJ_DIR}/qlog/idx_bam");

# step 7: featurecount: 
    %str_to_replace=();

    open IN, "<$FC_TEMPLATE";
    $content = do { local $/; <IN> };
    close IN;
    $str_to_replace{COMMON_SCRIPT_DIR} = $common_script_dir;
    $str_to_replace{PROJECT_DIR} = $proj_list{$proj}{PROJ_DIR};
    $str_to_replace{PROJECTNAME} = $proj;
    for my $k (keys %str_to_replace) {
        $content =~ s/\{$k\}/$str_to_replace{$k}/g;
    }

    open OUT, ">$FC_QSUB";
    print OUT $content;
    close OUT;

    system("mkdir -p $proj_list{$proj}{PROJ_DIR}/qlog/fc");

# step 8: deseq2: 
    %str_to_replace=();

    open IN, "<$DS2_TEMPLATE";
    $content = do { local $/; <IN> };
    close IN;
    
    open IN2, "<$DS2_R_TEMPLATE";
    my $content2 = do {local $/; <IN2>};
    close IN2;

    $str_to_replace{COMMON_SCRIPT_DIR} = $common_script_dir;
    $str_to_replace{PROJECT_DIR} = $proj_list{$proj}{PROJ_DIR};
    $str_to_replace{PROJECTNAME} = $proj;
    $str_to_replace{CONFIG_ROOT} = $proj_list{$proj}{CONFIG_ROOT};    
    $str_to_replace{DESEQ2_CONF_DIR} = $proj_list{$proj}{DESEQ2_CONF_DIR};
    $str_to_replace{DESEQ2_CONF} = $proj_list{$proj}{DESEQ2_CONF};
    $str_to_replace{DESEQ2_CONF_NUM} = $proj_list{$proj}{DESEQ2_CONF_NUM};
    $str_to_replace{DESEQ2_CONF_LIST} = $proj_list{$proj}{DESEQ2_CONF_LIST};

    for my $k (keys %str_to_replace) {
    	$content =~ s/\{$k\}/$str_to_replace{$k}/g;
    	$content2 =~ s/\{$k\}/$str_to_replace{$k}/g;
    }

    open OUT, ">$DS2_QSUB";
    print OUT $content;
    close OUT;

    open OUT2, ">$DS2_R";
    print OUT2 $content2;
    close OUT2;

    system("mkdir -p $proj_list{$proj}{PROJ_DIR}/qlog/deseq2");

# step 9: GSEA: 
    %str_to_replace=();

    open IN, "<$GSEA_TEMPLATE";
    $content = do { local $/; <IN> };
    close IN;
    open IN2, "<$GSEA_R_TEMPLATE";
    my $content2 = do { local $/; <IN2> };
    close IN2;
    
    $str_to_replace{COMMON_SCRIPT_DIR} = $common_script_dir;
    $str_to_replace{PROJECT_DIR} = $proj_list{$proj}{PROJ_DIR};
    $str_to_replace{PROJECTNAME} = $proj;
    $str_to_replace{DESEQ2_CONF_NUM} = $proj_list{$proj}{DESEQ2_CONF_NUM};
    $str_to_replace{DESEQ2_CONF_LIST} = $proj_list{$proj}{DESEQ2_CONF_LIST};
    $str_to_replace{SGE_TOTAL_TASK} = $proj_list{$proj}{SGE_TOTAL_TASK};
    for my $k (keys %str_to_replace) {
    	$content =~ s/\{$k\}/$str_to_replace{$k}/g;
	    $content2 =~ s/\{$k\}/$str_to_replace{$k}/g;
    }

    open OUT, ">$GSEA_QSUB";
    print OUT $content;
    close OUT;

    open OUT2, ">$GSEA_R";
    print OUT2 $content2;
    close OUT2;

    system("mkdir -p $proj_list{$proj}{PROJ_DIR}/qlog/gsea");
    
    
# now generate customized pipeline script
    %str_to_replace=();

    open IN, "<$PIPE_TEMPLATE";
    $content = do { local $/; <IN> };
    close IN;
    $str_to_replace{PROJECTNAME} = $proj;
    for my $k (keys %str_to_replace) {
	$content =~ s/\{$k\}/$str_to_replace{$k}/g;
    }

    open OUT, ">$PIPE_SCRIPT";
    print OUT $content;
    close OUT;

    # set the file permission to group readable:
    system("chmod -R 755 $proj_list{$proj}{PROJ_DIR}");
	
    print "The project $proj is created in $proj_list{$proj}{PROJ_DIR}\n";

} # end of each project

print "DONE!\n";


