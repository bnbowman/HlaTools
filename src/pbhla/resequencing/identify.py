import logging, os

log = logging.getLogger(__name__)

def identify_resequencing_data( summary, cns_fofn, read_fofn ):
    log.info("Identifying files for resequencing")
    selected = list( _parse_summary_file( summary ))
    cns_files = list( _parse_clusense_fofn( cns_fofn ))
    read_files = list( _parse_clusense_fofn( read_fofn ))
    resequencing_data = list( _select_resequencing_files( selected, cns_files, read_files ))
    log.info('Finished the idenifying the requisite files')
    return resequencing_data

def _parse_summary_file( summary ):
    log.info("Identifying files for resequencing based on selected contigs")
    with open( summary, 'r' ) as handle:
        handle.next() # Skip the header
        for line in handle:
            contig = line.split()[1]
            if contig == '-':
                continue
            yield contig

def _parse_clusense_fofn( fofn_file ):
    with open( fofn_file, 'r') as handle:
        for line in handle:
            filename = line.strip()
            if filename:
                yield filename

def _select_resequencing_files( selected, cns_files, read_files ):
    log.info('Finished the idenifying the requisite files')
    for cns_fn, read_fn in zip(cns_files, read_files):
        base_name = os.path.basename( read_fn )
        root_name = base_name.split('.')[0]
        if root_name in selected:
            yield (root_name, cns_fn, read_fn) 
