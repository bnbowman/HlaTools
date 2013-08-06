import os, logging

log = logging.getLogger()

class ContigInfo(object):
    def __init__(self, locus, contig, length, count, hit, pctid):
        self.locus = locus
        self.contig = contig
        self.length = int(length) if length.isdigit() else length
        self.count = int(count) if count.isdigit() else count
        self.hit = hit
        self.pctid = float(pctid)

    def __str__(self):
        return '\t'.join( [self.locus,
                           self.contig,
                           str(self.length),
                           str(self.count),
                           self.hit,
                           str(self.pctid) ])

def meta_summarize_contigs( contig_files, metadata, output_file, excluded=[] ):
    with open(output_file, 'w') as output:
        # Write metadata
        for key, value in metadata.iteritems():
            print >> output, "# %s=%s" % (key, value)
        # Write the column headers
        print >> output, "Locus\tContig\tLength\tCount\tHit\tPctId"
        # Finally parse and write the data
        for filepath in sorted(contig_files):
            basename = os.path.basename( filepath )
            locus = filepath.split('_')[-2]
            if locus in excluded:
                continue
            log.info( 'Reading the summary of %s from "%s"' % (locus, basename))
            summaries = list( _parse_contig_summaries( filepath, locus ))
            first, second = _select_contig_summaries( summaries, locus )
            print >> output, first
            print >> output, second

def _parse_contig_summaries( filepath, locus ):
    with open(filepath, 'r') as handle:
        handle.next() # Skip header
        for line in handle:
            info = [locus] + line.strip().split()
            yield ContigInfo( *info )

def _select_contig_summaries( summaries, locus ): 
    dummy = ContigInfo(locus, '-\t\t\t', '-', '-', '-', '0.0')
    if len( summaries ) == 0:
        return [dummy, dummy]
    elif len( summaries ) == 1:
        return [summaries[0], dummy]
    else:
        for summary in summaries:
            if summary.hit != summaries[0].hit:
                return [ summaries[0], summary ]
        return [ summaries[0], dummy ]
