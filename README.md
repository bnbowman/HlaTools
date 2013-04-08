# HlaTools #

HlaTools is a python package of tools for analyzing and phasing large
amplicons of genomic sequence.  Primarily this has been developed and
tested on the HLA genes from the human MHC region, but the same approach
should work for most complex regions of the human genome where phasing
is important.

## Requirements ##

The central pipeline in HlaTools is built upon HGAP, which in turn requires
a large portion of the PacBio Secondary Analysis Suite - Blasr, Quiver, etc.  

Currently the tools are tightly bound not just to HGAP, but to a modified 
version of HGAP installed locally on the PacBio servers.  Though currently
inconvenient, the hope is that we will be dis-entagled and given more sane 
and manageable dependencies in the near future.
