# HlaTools #

HlaTools is a python package of tools for analyzing and phasing large
amplicons of genomic sequence.  Primarily this has been developed and
tested on the HLA genes from the human MHC region, but the same approach
should work for most complex regions of the human genome where phasing
is important.

## Requirements ##

The central pipeline in HlaTools is built upon HBAR, formerly known as
HGAP, which is available online at:
http://github.com/PacificBiosciences/HBAR-DTK.git
HBAR-DTK itself as a fairly intensive installation procedure, detailed
in the README on the

In addition, the HlaTools pipeline also wraps two commandline utilites:
Blasr for the alignment of reads and Cd-Hit-Est for the reduction of
redundancy:
https://github.com/PacificBiosciences/blasr.git
http://weizhong-lab.ucsd.edu/cd-hit/

Currently the pipeline also has an option, still in development, to
resequence the resulting contigs with the PacBio Quiver consensus
algorithm.  If you wish to use this option, you will also need the
full PacBio SMRT Analysis installation available here:
http://pacbiodevnet.com
