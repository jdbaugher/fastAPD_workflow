fastAPD_workflow
================

### Nucleotide diversity (average pairwise difference) workflow

This workflow makes use of the Bio::fastAPD Perl module for fast calculations
of nucleotide diversity (average pairwise difference (APD)) from 
multiple sequence alignments. The workflow accompanies the following peer-reviewed 
journal article:

Julian TR, Baugher JD, Rippinger CM, Pinekenstein R, Kolawole AO, Mehoke TS, Wobus CE, 
Feldman AB, Pineda FJ, Schwab KJ. Murine norovirus (MNV-1) exposure in vitro to the the purine
nucleoside analog Ribavirin increases quasispecies diversity. Virus research. 2016 Jan 4;211:165-173.

## Requirements
This workflow has been tested on MACOSX and linux operating systems
using recent versions of Perl, R, and X11.

#### Perl modules:
    Bio::fastAPD v1.10.0 or higher
#### R libraries:
    ggplot2

### Data for Analysis
A subdirectory named <i>data</i> may be created in the <i>fastAPD_workflow</i> directory to house
the data for analysis.

## Usage

To run the analysis - 
From within a shell session, navigate to the <i>fastAPD_workflow</i> directory and type:

    bash nucleotide_diversity.sh ./data

The analysis consists of four steps. The user will be notified when the script 
has finished. If no errors are reported, the analysis is complete.


## Authors

Joseph D. Baugher, Ph.D. and Fernando J. Pineda, Ph.D.<br>
Copyright (c) 2014,2015 Joseph D. Baugher, Ph.D.

## Maintainer

Joseph D. Baugher, Ph.D., joebaugher(at)hotmail.com
