import re
import csv

from pbcore.io.FastaIO import FastaReader

HEADER = ["Group", "Reads",       "FpReads",
                   "FlReads",     "Bp",         
                   "FpBp",        "FlBp",
                   "AvgLength",   "FpAvgLength",
                   "FlAvgLength", "FpFraction", 
                   "FlFraction",  "ReadFraction"]

MARGIN = 15

class SequencingStats( object ):
    """
    A class for tabulating subread statistics for a single HLA Locus
    """

    def __init__(self, locus, start=None, end=None):
        self.locus = locus
        self.start = int(start) + MARGIN if start else None
        self.end = int(end) - MARGIN if end else None
        self.all_subreads = 0
        self.subreads = 0
        self.fullpass_subreads = 0
        self.full_length_subreads = 0
        self.bases = 0
        self.fullpass_bases = 0
        self.full_length_bases = 0

    @property
    def average_length(self):
        if self.subreads:
            return self.bases/self.subreads
        return 0

    @property
    def average_fullpass_length(self):
        if self.fullpass_subreads:
            return self.fullpass_bases/self.fullpass_subreads
        return 0

    @property
    def average_full_length(self):
        if self.full_length_subreads:
            return self.full_length_bases/self.full_length_subreads
        return 0

    @property
    def fullpass_fraction(self):
        if self.subreads:
            return round(self.fullpass_subreads/float(self.subreads), 3)
        return 0

    @property
    def full_length_fraction(self):
        if self.subreads:
            return round(self.full_length_subreads/float(self.subreads), 3)
        return 0

    @property
    def subread_fraction(self):
        if self.subreads:
            return round(self.subreads/float(self.all_subreads), 3)
        return 0

    def alignment_overlaps_amplicon(self, alignment):
        if self.start is None or self.end is None:
            return False
        if self.start >= alignment.pos and self.end <= alignment.aend:
            return True
        return False

    def add_aligned_read(self, alignment):
        self.subreads += 1
        self.bases += int(alignment.tlen)
        if alignment.qname.startswith('fp'):
            self.fullpass_subreads += 1
            self.fullpass_bases += alignment.tlen
        if self.alignment_overlaps_amplicon( alignment ):
            self.full_length_subreads += 1
            self.full_length_bases += alignment.tlen

    def tabulate(self):
        return [self.locus, self.subreads,             self.fullpass_subreads,
                            self.full_length_subreads, self.bases,             
                            self.fullpass_bases,       self.full_length_bases,
                            self.average_length,       self.average_fullpass_length,
                            self.average_full_length,  self.fullpass_fraction,    
                            self.full_length_fraction, self.subread_fraction]

class SubreadStats( object ):
    """
    A class for tabulating subread statistics for the HLA Pipeline
    """

    def __init__(self, ref_fasta, locus_dict, amplicon_csv=None):
        self.loci = {}
        self.references = {}
        self.amplicons = {}
        self.amplicon_locations = {}
        self.locus_dict = locus_dict
        self.all_stats = SequencingStats( 'All' )
        if amplicon_csv is None:
            for record in FastaReader(ref_fasta):
                self.add_stats( record.name )
        else:
            with open(amplicon_csv, 'r') as handle:
                handle.readline() # Skip header line
                for row in csv.reader(handle):
                    ref, locus, amplicon, start, end = row
                    self.add_stats( ref, start, end )

    def add_stats(self, name, start=None, end=None, amplicon=None):
        # Add a high-level, locus-wide object
        locus = self.locus_dict[name]
        if locus not in self.loci:
            self.loci[locus] = SequencingStats( locus )
        # Add a mid-level, gene-specific object
        if name not in self.references:
            self.references[name] = SequencingStats( name )
        # Add a low-level, amplicon-specific object
        amplicon_name = "{0}_{1}".format(name, amplicon)
        if amplicon_name not in self.amplicons:
            self.amplicons[amplicon_name] = SequencingStats( amplicon_name, start, end )
            try:
                self.amplicon_locations[name].append( (amplicon_name, start, end) )
            except:
                self.amplicon_locations[name] = [( amplicon_name, start, end )]

    def find_best_amplicon(self, alignment):
        if len(self.amplicon_locations[alignment.rname]) == 1:
            return self.amplicon_locations[alignment.rname][0][0]
        results = {}
        for amp, amp_start, amp_end in self.amplicon_locations[alignment.rname]:
            overlap = calculate_overlap( alignment.pos, alignment.aend, amp_start, amp_end )
            results[amp] = overlap
        return sorted(results)[0][0]

    def add_aligned_read( self, locus, alignment ):
        self.all_stats.add_aligned_read( alignment )
        self.loci[locus].add_aligned_read( alignment )
        self.references[alignment.rname].add_aligned_read( alignment )
        amplicon = self.find_best_amplicon( alignment )
        self.amplicons[amplicon].add_aligned_read( alignment )
        self.increment_readtotals()

    def increment_readtotals(self):
        self.all_stats.all_subreads += 1
        for locus, stats in self.loci.iteritems():
            stats.all_subreads += 1
        for reference, stats in self.references.iteritems():
            stats.all_subreads += 1
        for amplicon, stats in self.amplicons.iteritems():
            stats.all_subreads += 1

    def write_summary(self, handle, data):
        writer = csv.writer( handle )
        writer.writerow( HEADER )
        for name, stats in data.iteritems():
            writer.writerow( stats.tabulate() )
        writer.writerow( self.all_stats.tabulate() )

    def write(self, output_file, group='locus'):
        with open(output_file, "w") as handle:
            if group == 'locus':
                self.write_summary(handle, self.loci )
            if group == 'reference':
                self.write_summary(handle, self.references )
            if group == 'amplicon':
                self.write_summary(handle, self.amplicons )
            
#
# Utilities
#
def calculate_overlap(start1, end1, start2, end2):
    if start1 >= end2:
        return 0
    elif end1 <= start2:
        return 0
    elif start1 <= start2 and end1 >= end2:
        return end2-start2
    elif start1 >= start2 and end1 <= end2:
        return end1-start1
    elif start1 <= start2 and end1 <= end2:
        return end1-start2
    elif start1 >= start2 and end1 >= end2:
        return end2-start1
