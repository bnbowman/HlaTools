
import logging
import md5
import os

from pbcore.io import FastaReader

log = logging.getLogger()

def fasta_movie_counts( fasta ):
    counts = {'all':0}
    for record in FastaReader( fasta ):
        movie = record.name.split('_')[0]
        counts['all'] += 1
        try:
            counts[movie] += 1
        except:
            counts[movie] = 1
    return counts

def fofn_file_counts( fasta_fofn ):
    counts = {'all':0}
    with open( fasta_fofn ) as handle:
        for line in handle:
            filename = line.strip()
            movie_name = ".".join( os.path.basename(filename).split(".")[:-1] )
            md5_name = md5.md5(movie_name).hexdigest()[:8]
            for record in FastaReader( filename ):
                movie = record.name.split('/')[0]
                counts['all'] += 1
                try:
                    counts[md5_name] += 1
                except:
                    counts[md5_name] = 1
    return counts

def fofn_naming_dict( fasta_fofn ):
    names = {}
    with open( fasta_fofn ) as handle:
        for line in handle:
            filename = line.strip()
            movie_name = ".".join( os.path.basename(filename).split(".")[:-1] )
            md5_name = md5.md5(movie_name).hexdigest()[:8]
            for record in FastaReader( filename ):
                name = record.name.split()[0]
                movie, hole, pos = name.split('/')[0:3]
                start = pos.split('_')[0]
                short_name = '%s_%s_%s' % (md5_name, hole, start)
                names[short_name] = name
    return names

def write_renaming_key(subread_fofn, renamed_subreads, output_file ):
    """
    Create a key for translating HBAR subread names to canonical PacBio names
    """
    # Compare the two files to make sure they're equivalent
    raw_counts = fofn_file_counts( subread_fofn )
    new_counts = fasta_movie_counts( renamed_subreads )
    try:
        assert raw_counts == new_counts
    except AssertionError:
        msg = 'The number of raw subreads (%s) does not ' % raw_count + \
              'match the number of renamed reads (%s)' % new_count
        log.info( msg )
        raise ValueError( msg )
    # Write out the pairs of names to file
    name_dict = fofn_naming_dict( subread_fofn )
    with open( output_file, 'w') as handle:
        for record in FastaReader(renamed_subreads):
            raw_name = record.name.split()[0]
            new_name = name_dict[raw_name]
            handle.write('%s\t%s\n' % (raw_name, new_name))
    return output_file

if __name__ == '__main__':
    import sys

    subread_fofn = sys.argv[1]
    renamed_subreads = sys.argv[2]
    output_file = sys.argv[3] if len(sys.argv) > 3 else sys.stdout

    write_renaming_key( subread_fofn, renamed_subreads, output_file )
