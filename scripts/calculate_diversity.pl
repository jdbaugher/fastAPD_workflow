#!/usr/bin/env perl -w
# calculate_diversity.pl
# Author: Joseph D. Baugher, <jbaughe2(at)jhmi.edu>
# Copyright (c) 2014 Joseph D. Baugher

use strict;
use Getopt::Long;
use File::chdir;
use File::Spec;

use Bio::fastAPD;

my $input_dir = "./";
my $glob      = "*";
my $prefix    = "";                 

my $the_options = GetOptions (
    "input_dir=s" => \$input_dir,
    "glob=s"      => \$glob,
    "prefix=s"    => \$prefix    
    );

print join("\t", qw/ID Lab Amplicon Group Passage Replicate  
                     fastAPD Std_Err N_Pos N_Reads/), "\n"; 
                    
{
    local $CWD = $input_dir;
    my @input_files = glob($glob);

    foreach my $ifile (@input_files) {

        my ($volume, $path, $filename) = File::Spec->splitpath($ifile);
        my $subject = join("_", $prefix, (split("_",$filename))[0]);
                        
        ################ Project specific annotation ##################
        
        my @info = qw/NA Stock_P0_Lab2_R1 Stock_P0_Lab2_R2 Stock_50RV_P0_Lab2_R1
                Stock_50RV_P0_Lab2_R2 Control_P5_Lab2_R1 Control_P5_Lab2_R2 
                Control_P5_Lab2_R3 Ribavirin_50RV_P5_Lab2_R1 
                Ribavirin_50RV_P5_Lab2_R2 Ribavirin_50RV_P5_Lab2_R3
                Control_P10_Lab2_R1 Control_P10_Lab2_R2 Control_P10_Lab2_R3 
                Ribavirin_50RV_P10_Lab2_R1 Ribavirin_50RV_P10_Lab2_R2
                Ribavirin_50RV_P10_Lab2_R3 Control_P15_Lab2_R1 Control_P15_Lab2_R2
                Control_P15_Lab2_R3 Ribavirin_50RV_P15_Lab2_R1
                Ribavirin_50RV_P15_Lab2_R2 Ribavirin_50RV_P15_Lab2_R3 
                Control_P20_Lab2_R1 Control_P20_Lab2_R2 Control_P20_Lab2_R3
                Ribavirin_50RV_P20_Lab2_R1 Ribavirin_50RV_P20_Lab2_R2 NA 
                Control_P5_Lab1_R1 Control_P5_Lab1_R2 Control_P5_Lab1_R3
                Ribavirin_100RV_P5_Lab1_R1 Ribavirin_100RV_P5_Lab1_R2 
                Ribavirin_100RV_P5_Lab1_R3 Control_P10_Lab1_R1
                Control_P10_Lab1_R2 Control_P10_Lab1_R3 
                Ribavirin_100RV_P10_Lab1_R1 Ribavirin_100RV_P10_Lab1_R2
                Ribavirin_100RV_P10_Lab1_R3 Control_P15_Lab1_R1 
                Control_P15_Lab1_R2 Control_P15_Lab1_R3
                Ribavirin_100RV_P15_Lab1_R1 Ribavirin_100RV_P15_Lab1_R2 
                Ribavirin_100RV_P15_Lab1_R3 Control_P20_Lab1_R1
                Control_P20_Lab1_R2 Control_P20_Lab1_R3 
                Ribavirin_100RV_P20_Lab1_R1 Ribavirin_100RV_P20_Lab1_R2
                Ribavirin_100RV_P20_Lab1_R3/; 

        # Capture the 2-digit ID number 
        # e.g. id = 14 for file name 
        #    APL000000014.78.v0p2_Amp3.genotypes.sig.all.expanded_haplotypes.fa
        $_ = (split("\\.", $ifile))[0];
        /([1-9]?[0-9])$/;
        my $id = $1;

        # Use the id as an index for the @info array which contains additional info
        # Store the laboratory
        my $lab;    
        if ($info[$id] =~ /Lab1/) {$lab = "1"}
        elsif ($info[$id] =~ /Lab2/) {$lab = "2"}
        else{$lab = "NA"}
        
        # Store the amplicon number
        $_ = $ifile;
        /Amp([2-4])/;
        my $amplicon = $1;
        
        # Store the passage number
         $_ = $info[$id];
        /P(\d+)/;
        my $passage = $1;
        
        # Store the group (Treatment or Control)
        my $group;
        if ($info[$id] =~ /Ribavirin/) {$group = "T"}
        elsif ($info[$id] =~ /Control|Stock/) {$group = "C"}
        else{$group = "NA"}            
            
        # Store the replicate number    
        $_ = $info[$id];
        /(\d)$/;
        my $replicate = $1;            
        
        ######## Mask positions with outlier mutation rates
        my $mask;
        if ($amplicon == 2) {$mask = (1 x 420)}
        elsif ($amplicon == 3) {
            # Exclude Amp3 positions 11 and 292
            $mask = join("", (1 x 10),0,(1 x (291-11)),0,(1 x (376-292)));
        }
        elsif ($amplicon == 4) {
            # Exclude Amp4 positions 326, 361, 393
            $mask = join("",(1 x 325),0,(1 x (360-326)),0,(1 x (392-361)),0,(1 x (397-393)));
        }

        #################### Perform Calculations #################

        ### Create an array of sequences from an MSA file
        my $sequences_ref = create_seq_array($ifile);   
	my $fastAPD_obj = Bio::fastAPD->new();
        $fastAPD_obj->initialize(seq_array_ref => $sequences_ref,
                                 alphabet      => 'dna'
                                 mask          => $mask);
                                 
        # Perform calculations                         
        my $apd = $fastAPD_obj->calculate_diversity(method  => 'fast_apd',
                                            	    compare => 'gap_base');                                      
        my $std_err = $fastAPD_obj->calculate_apd_std_err(method  => 'fast_apd',
                                                  	  compare => 'gap_base');                          
        my $num_reads     = $fastAPD_obj->n_reads;                                      
        my $num_positions = $fastAPD_obj->width; 
          
        print join("\t", $id, $lab, $amplicon, $group, $passage, $replicate,
            $apd, $std_err, $num_positions, $num_reads), "\n";                
    }
}

####################### Subroutines ###########################

sub create_seq_array {
    my $file = shift;
    my @sequences;
                                                                       
    # Read in aligned sequence file (fasta format for this test)
    open(my $input_fh, '<', $file);
    chomp(my @fasta_lines = <$input_fh>);
    close $input_fh;
     
    # Create an array of aligned sequences
    my $curr_seq;
    foreach my $line (@fasta_lines) {
        if (substr($line, 0, 1) eq ">") { 
            if ($curr_seq) { push(@sequences, $curr_seq) }
            $curr_seq = ();
        }
        else { $curr_seq .= $line }
    }
    if ($curr_seq) { push(@sequences, $curr_seq) }
 
    return(\@sequences);
}

__END__

