# HlaTools #

HlaTools is a python package of tools for analyzing and phasing large
amplicons of genomic sequence.  Primarily this has been developed and
tested on the HLA genes from the human MHC region, but the same approach
should work for most complex regions of the human genome where phasing
is important.

## Installation ##

For the primary tools provided by the HlaTools package, we require 5 other
python packages as dependencies:
`Numpy:        http://www.numpy.org/`
`PBcore:       http://github/PacificBiosciences/pbcore.git`
`PypeFlow:     http://github/cschin/pypeFLOW.git`
`PBdagcon:     http://github/PacificBiosciences/pbdagcon.git`
`PhasingTools: http://github/bnbowman/PhasingTools.git`

All 5 packages can be quickly and easily installed with pip:
`$ pip install numpy             # For NumPy, which is in PyPi`
`$ pip install git+<GitHub URL>  # For all other packages`

## Optional Requirements ##

The other major script included in HlaTools is a pipeline for de-convoluting
complex HLA amplicon samples with a mixture of Class I and Class II genes
or amplicons with large discrepancies in length.  This pipeline is built
upon HBAR, formerly known as HGAP, and therefore also requires the HBAR
development package in order to use.  The HBAR-DTK itself has a fairly
intensive installation procedure, detailed in the README distributed with
the tool from the HBAR-DTK github page:
`http://github.com/PacificBiosciences/HBAR-DTK.git`

In addition, this pipeline also wraps two simpler commandline utilites:
Blasr for the alignment of reads and Cd-Hit-Est for the reduction of
redundancy:
`https://github.com/PacificBiosciences/blasr.git`
`http://weizhong-lab.ucsd.edu/cd-hit/`