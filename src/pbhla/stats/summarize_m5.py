#! /usr/bin/env python

import sys

from pbhla.io.BlasrIO import BlasrReader

def summarize_m5( m5_file, output_file ):
    with open( output_file, 'w' ) as handle:
        handle.write("Query\tRef\tMatch\tMis\tIns\tDel\tErr\n")
        for entry in BlasrReader( m5_file ):
            # Count the bases of each type
            nmat = int(entry.nmat)
            nmis = int(entry.nmis)
            nins = int(entry.nins)
            ndel = int(entry.ndel)
            nerr = nmis + nins + ndel
            ntot = nmat + nerr
            # Calculate the ratios of each type
            rmat = round(nmat / float(ntot), 4)
            rmis = round(nmis / float(ntot), 4)
            rins = round(nins / float(ntot), 4)
            rdel = round(ndel / float(ntot), 4)
            rerr = round(nerr / float(ntot), 4)
            # Write the output to file
            line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ( entry.qname,
                                                      entry.tname, 
                                                      rmat,
                                                      rmis,
                                                      rins,
                                                      rdel,
                                                      rerr )
            handle.write( line )

if __name__ == '__main__':
    m5_file = sys.argv[1]
    output_file = sys.argv[2]

    summarize_m5( m5_file, output_file )
