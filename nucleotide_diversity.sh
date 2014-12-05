#!/bin/bash
#######################
# nucleotide_diversity.sh
# Author: Joseph D. Baugher, jbaughe2(at)jhmi.edu
# Copyright © 2014 Joseph D. Baugher
#######################
#$ -N nucleotide_diversity
#$ -cwd
#$ -S /bin/bash 
#######################

# If starting with haplotypes, they must expanded based on the number of reads 
# represented by each haplotype. Otherwise, skip this step.
echo "Running Step 1 - haplotype expansion..."
perl scripts/expand_haplotypes.pl --input_dir $1 --glob *genotypes.sig.txt

# The nucleotide diversity values are estimated
echo "Running Step 2 - nucleotide diversity calculation..."
perl scripts/calculate_diversity.pl --input_dir $1 --glob *expanded_haplotypes.fa > diversity_results.txt

# The amplicon results are plotted
echo "Running Step 3 - plotting results..."
Rscript scripts/plot_Ribavirin_diversity.Rscript

# Calculate stats
echo "Running Step 4 - stats..."
Rscript scripts/diversity_stats.Rscript
echo "Analysis finished."
