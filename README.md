fastAPD_workflow
================

### Nucleotide diversity (average pairwise difference) workflow

This workflow makes use of the Bio::fastAPD Perl module for fast calculations
of nucleotide diversity (average pairwise difference (APD)) from 
multiple sequence alignments. The workflow accompanies an upcoming peer-reviewed 
journal article. Additional details will be provided as they become available.

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
