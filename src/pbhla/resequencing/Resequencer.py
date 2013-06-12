#! /usr/bin/env python

import os
import logging

from pbhla.utils import create_directory

log = logging.getLogger()

def resequence( smrt_analysis, raw_data_file, reference_file, output_dir ):
    ### TODO: maybe adjust compare seq options so that euqal mapping are not assigned randomly...
    ### we resequence once, elminate redundant clusters, then resequence the remaining clusters
    ### to avoid splitting our CLRs over clusters that respresent the same template 
    log.info('Resequencing "{0}"')
    create_directory( output_dir )

    ### create a .cmp.h5
    reference_alignment = os.path.join( output_dir, "reference.cmp.h5")
    if os.path.isfile( reference_alignment ):
        log.info('Found existing reference alignment "%s"' % reference_alignment)
        log.info('Skipping reference alignment step...')
    else:
        log.info('No existing reference alignment found, creating...')
        smrt_analysis.compare_sequences( raw_data_file, 
                                         reference_file,
                                         reference_alignment )
        smrt_analysis.load_pulses( raw_data_file, reference_alignment )
        smrt_analysis.sort_cmpH5( reference_alignment )
        log.info('No existing reference alignment found, creating...')
    
    ### run quiver 
    consensus_fasta = os.path.join(output_dir, "consensus.fasta")
    consensus_fastq = os.path.join(output_dir, "consensus.fastq")
    if os.path.isfile( consensus_fasta ):
        log.info('"Found existing Quiver Results "{0}"'.format(consensus_fasta))
        log.info('Skipping Quiver resequencing step...')
    else:
        log.info('No existing quiver output found, creating...')
        smrt_analysis.quiver( reference_alignment,
                              reference_file, 
                              consensus_fasta, 
                              consensus_fastq)
        log.info('Finished Quiver resequencing')

    return consensus_fasta
    
"""
    ### separate out resequenced HLA sequences from the nonspecific ones, by sequence name, store in two different files
    ### TODO: read fastq output, parse consnesus QV, store somewhere the annotator can get at it
    target_seq_names=[]
    for record in FastaReader(self.hap_con_fasta):
        target_seq_names.append( record.name )
    with FastaWriter(self.args.proj+"/reseq/resequenced_hap_con.fasta") as on_target: 
        with FastaWriter(self.args.proj+"/reseq/resequenced_nonspecific.fasta") as off_target:
            for record in FastaReader(reseq_output+".fasta"):
                if record.name.split("|")[0] in target_seq_names:
                    on_target.writeRecord( record )
                else:
                    off_target.writeRecord( record )
    with FastqWriter(self.args.proj+"/reseq/resequenced_hap_con.fastq") as on_target:
        with FastqWriter(self.args.proj+"/reseq/resequenced_nonspecific.fastq") as off_target:
            for record in FastqReader(reseq_output+".fastq"):
                if record.name.split("|")[0] in target_seq_names:
                    on_target.writeRecord( record )
                else:   
                    off_target.writeRecord( record )
"""
