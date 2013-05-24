#!/home/UNIXHOME/jquinn/HGAP_env/bin/python

import logging

from pbhla.io.SamIO import SamReader

BIN_SIZE = 50

class AmpliconFinder( object ):
    """
    A class for tabulating subread statistics for a single HLA Locus
    """

    def __init__(self, sam_file, output_file, locus_dict):
        self.sam_file = sam_file
        self.output_file = output_file
        self.locus_dict = locus_dict
        self.initialize_logger()
        self.run()

    def initialize_logger(self):
        self.log = logging.getLogger(__name__)
        self.log.info("Initializing AmpliconFinder with:")
        self.log.info("\tInput Data: {0}".format(self.sam_file))
        self.log.info("\tOutput Data: {0}".format(self.output_file))

    def run(self):
        self.parse_positions()
        self.identify_amplicons()
        self.write_output()

    def parse_positions(self):
        self.starts = {}
        self.ends = {}
        for entry in SamReader( self.sam_file ):
            try:
                self.starts[entry.rname].append( entry.pos )
            except:
                self.starts[entry.rname] = [ entry.pos ]
            try:
                self.ends[entry.rname].append( entry.aend )
            except:
                self.ends[entry.rname] = [ entry.aend ]

    def identify_amplicons(self):
        self.amplicons = []
        for ref in self.starts:
            locus = self.locus_dict[ref]
            if locus in ['DQB']:
                count = 2
            else:
                count = 1
            start_peaks = self.find_peaks(self.starts[ref], count)
            end_peaks = self.find_peaks(self.ends[ref], count)
            pos_pairs = zip(start_peaks, end_peaks)
            for i, pair in enumerate(pos_pairs):
                start, end = pair
                amplicon = (ref, locus, str(i+1), str(start), str(end))
                self.amplicons.append( amplicon )
            
    def find_peaks(self, positions, count):
        binned_pos = [p/BIN_SIZE*BIN_SIZE for p in positions]
        min_bin = -1 * BIN_SIZE
        max_bin = max(positions) + 2 * BIN_SIZE
        bins = range(min_bin, max_bin, BIN_SIZE)
        counts = [binned_pos.count(b) for b in bins]
        scores = [counts[c]-max(counts[c-1], counts[c+1])
                         for c in range(1,len(counts)-1)]
        sort_scores = sorted( scores, reverse=True )
        peak_bin_pos = [scores.index(s)+1 for s in sort_scores[:count]]
        peak_bins = [bins[p] for p in peak_bin_pos]
        peak_ranges = [(p-BIN_SIZE, p+2*BIN_SIZE) for p in peak_bins]
        peak_pos = [subset_positions(positions, start, end) for start, end in peak_ranges]
        peak_loc = [median(p) for p in peak_pos]
        return sorted(peak_loc)

    def write_output(self):
        with open(self.output_file, "w") as handle:
            print >> handle, "Reference,Locus,Amplicon,Start,End"
            print self.amplicons
            for amplicon in sorted(self.amplicons, key=lambda x: x[1]):
                print >> handle, ','.join( amplicon )

#
# Utilities
#
def subset_positions(positions, start, end):
    return [p for p in positions
                    if p >= start
                    if p <= end]

def median(positions):
    positions = sorted(positions)
    if len(positions) % 2:
        middle = len(positions)/2+1
    else:
        middle = len(positions)/2
    return positions[middle]
