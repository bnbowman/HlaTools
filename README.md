# HlaTools #

HlaTools is a python package of tools for analyzing and phasing large
amplicons of genomic sequence.  Primarily this has been developed and
tested on the HLA genes from the human MHC region, but the same approach
should work for most complex regions of the human genome where phasing
is important.

## Installation ##

For the primary tools provided by the HlaTools package, we require 5 other
python packages as dependencies:<br>
`Numpy:        http://www.numpy.org/`<br>
`PBcore:       http://github/PacificBiosciences/pbcore.git`<br>
`PypeFlow:     http://github/cschin/pypeFLOW.git`<br>
`PBdagcon:     http://github/PacificBiosciences/pbdagcon.git`<br>
`PhasingTools: http://github/bnbowman/PhasingTools.git`<br>

All 5 packages can be quickly and easily installed with pip:<br>
`$ pip install numpy             # For NumPy, which is in PyPi`<br>
`$ pip install git+<GitHub URL>  # For all other packages`

## Optional Requirements ##

The other major script included in HlaTools is a pipeline for de-convoluting
complex HLA amplicon samples with a mixture of Class I and Class II genes
or amplicons with large discrepancies in length.  This pipeline is built
upon HBAR, formerly known as HGAP, and therefore also requires the HBAR
development package in order to use.  The HBAR-DTK itself has a fairly
intensive installation procedure, detailed in the README distributed with
the tool from the HBAR-DTK github page:<br>
`http://github.com/PacificBiosciences/HBAR-DTK.git`

In addition, this pipeline also wraps two simpler commandline utilites:
Blasr for the alignment of reads and Cd-Hit-Est for the reduction of
redundancy:<br>
`https://github.com/PacificBiosciences/blasr.git`<br>
`http://weizhong-lab.ucsd.edu/cd-hit/`