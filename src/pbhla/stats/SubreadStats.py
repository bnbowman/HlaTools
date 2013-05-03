#!/home/UNIXHOME/jquinn/HGAP_env/bin/python
import re
import csv

from pbhla.io.FastaIO import FastaReader

HEADER = ["Locus", "Reads",     "FpReads",
                   "Bp",        "FpBp",
                   "AvgLength", "FpAvgLength",
                   "FpFraction", "ReadFraction"]

class LocusStats( object ):
    """
    A class for tabulating subread statistics for a single HLA Locus
    """

    def __init__(self, locus ):
        self.locus = locus
        self._all_subreads = 0
        self._subreads = 0
        self._fullpass_subreads = 0
        self._bases = 0
        self._fullpass_bases = 0

    @property
    def subreads(self):
        return self._subreads

    @property
    def fullpass_subreads(self):
        return self._fullpass_subreads

    @property
    def bases(self):
        return self._bases

    @property
    def fullpass_bases(self):
        return self._fullpass_bases

    @property
    def average_length(self):
        return self._bases/self._subreads

    @property
    def average_fullpass_length(self):
        return self._fullpass_bases/self._fullpass_subreads

    @property
    def fullpass_fraction(self):
        return round(self._fullpass_subreads/float(self._subreads), 3)

    @property
    def subread_fraction(self):
        return round(self._subreads/float(self._all_subreads), 3)

    def increment_readtotal(self):
        self._all_subreads += 1

    def add_aligned_read(self, alignment):
        self._subreads += 1
        self._bases += int(alignment.tlen)
        if alignment.qname.startswith('fp'):
            self._fullpass_subreads += 1
            self._fullpass_bases += int(alignment.tlen)

    def tabulate(self):
        return [self.locus, self.subreads,          self.fullpass_subreads, 
                            self.bases,             self.fullpass_bases,
                            self.average_length,    self.average_fullpass_length,
                            self.fullpass_fraction, self.subread_fraction]

class SubreadStats( object ):
    """
    A class for tabulating subread statistics for the HLA Pipeline
    """

    def __init__(self, ref_fasta, locus_dict):
        self.loci = {}
        for record in FastaReader(ref_fasta):
            locus = locus_dict[record.name]
            self.loci[locus] = LocusStats( locus )
        self.loci['All'] = LocusStats( 'All' )

    def add_aligned_read( self, locus, alignment):
        self.loci[locus].add_aligned_read( alignment )
        self.loci['All'].add_aligned_read( alignment )
        self.increment_readtotals()

    def increment_readtotals(self):
        for locus, stats in self.loci.iteritems():
            stats.increment_readtotal()

    def tabulate_locus( self, locus ):
        return self.loci[locus].tabulate()

    def write(self, output_file):
        with open(output_file, "w") as handle:
            writer = csv.writer( handle )
            writer.writerow( HEADER )
            for name, stats in self.loci.iteritems():
                if name == 'All':
                    continue
                writer.writerow(self.tabulate_locus( name )) 
            writer.writerow(self.tabulate_locus( "All" ))
