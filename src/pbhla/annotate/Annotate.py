
import os, sys, argparse, logging

from collections import namedtuple

from pbcore.io.FastaIO import FastaReader
from pbcore.io.FastqIO import FastqReader
from pbcore.io.GffIO import GffWriter

from pbhla.io.GffIO import ( create_annotation,
                             create_var_annotation )
from pbhla.utils import create_directory
from pbhla.fasta.utils import write_fasta

log = logging.getLogger()

info = namedtuple('info', 'canonical_pos, feature, codon')

def annotate(self, file_list, output_dir):
    log.info("Initializing Annotation Process")
    create_directory( file_list )

    locus_dict={}
    with open(self.args.proj+"/phased/phasr_output_seqs.txt", "r") as f:
        for line in f:
            seq_name, locus = line.strip().split()
            self.locus_dict[seq_name] = locus
    log.info("Found %s sequences to annotate" % len(locus_dict))

    MSA_fn_dict={}; 
    MSA_cDNA_fn_dict={}
    MSA_info_fn_dict={}; 
    MSA_cDNA_info_fn_dict={}

    tmp_fn = os.path.join( output_dir, 'tmp.fasta' )
    with open(self.args.MSA, "r") as f:
        for line in f:
            genomic, nucleotide, locus = line.split()
            MSA_fn_dict[locus] = genomic + ".afa"
            MSA_info_fn_dict[locus] = genomic + ".info"
            MSA_cDNA_info_fn_dict[locus] = nucleotide + ".info"
            MSA_cDNA_fn_dict[locus] = nucleotide + ".afa"

    gDNA_var_path = os.path.join(output_dir, 'best_gDNA_ref_match_variants.gff')
    gDNA_var_writer = GffWriter( gDNA_var_path )
    gDNA_var_writer.writeMetaData('pacbio-variant-version', '1.4')

    cDNA_var_path = os.path.join(output_dir, 'best_cDNA_ref_match_variants.gff')
    cDNA_var_writer = GffWriter( cDNA_var_path )
    cDNA_var_writer.writeMetaData('pacbio-variant-version', '1.4')

    gDNA_annot_path = os.path.join(output_dir, "gDNA.gff")
    gDNA_annot_writer = GffWriter( gDNA_annot_path )

    cDNA_annot_path = os.path.join(output_dir, "cDNA.gff")
    cDNA_annot_writer = GffWriter( cDNA_annot_path )

    feature_path = os.path.join(output_dir, "features.gff")
    feature_writer = GffWriter( feature_path )

    gDNA_calls={}
    cDNA_calls={}

    gDNA_allele_path = os.path.join(output_dir, 'gDNA_allele_calls.txt')
    if os.path.isfile( gDNA_allele_path ):
        os.remove( gDNA_allele_path )

    cDNA_allele_path = os.path.join(output_dir, 'cDNA_allele_calls.txt')
    if os.path.isfile( cDNA_allele_path ):
        os.remove( cDNA_allele_path )

    hap_con_fasta_path = os.path.join(output_dir, 'resequenced_hap_con_cDNA.fasta')
    if os.path.isfile( hap_con_fasta_path ):
        os.remove( hap_con_fasta_path )

    hap_con_fastq_path = os.path.join(output_dir, 'resequenced_hap_con_cDNA.fastq')
    for r in FastqReader( hap_con_fastq_path ):
        write_fasta([r], tmp_fn, "w") 
        name = r.name.split("|")[0]
        locus = locus_dict[name]
        log.info('Processing "%s"  from locus "%s"' % (name, locus))

        ### read in profile features
        MSA_info_fn = MSA_info_fn_dict[locus]
        MSA_info = {}
        with open(MSA_info_fn, "r") as f2:
            for line in f2:
                pos1, pos2, region = line.strip().split()
                MSA_info[int(pos1)] = info._make( ( pos2, region, 0 ) )

        ### read profile fn
        MSA_fn = MSA_fn_dict[locus]

        ### align sequence to profile for this locus
        output = os.path.join(output_dir, name + ".afa")
        if not os.path.isfile(output):
            try:
                muscle_output = check_output(". /mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh; \
                    muscle -profile -in1 %s -in2 %s -out %s 2> /dev/null" % (MSA_fn, tmp_fn, output),
                    executable='/bin/bash', shell=True)
            except:
                log.info('MSA failed for "%s"' % name )
                os.remove( tmp_fn )
                continue

        ### read out a reference sequence from the old profile
        letters_max = 0
        for r2 in FastaReader(MSA_fn):
            num_letters = ( r2.sequence.count('A') + r2.sequence.count('G') + r2.sequence.count('C') + r2.sequence.count('T') )
            if num_letters > letters_max:
                letters_max = num_letters
        for r2 in FastaReader(MSA_fn):
            num_letters = ( r2.sequence.count('A') + r2.sequence.count('G') + r2.sequence.count('C') + r2.sequence.count('T') )
            if num_letters == letters_max:
                comparison_seq = r2.sequence
                comparison_seq_name = r2.name
                break

        ### read out the same reference sequence from the NEW profile
        ### and read out our consensus sequence
        for r2 in FastaReader(output):
            if r2.name == comparison_seq_name:
                new_comparison_seq = r2.sequence
                new_comparison_seq_name = r2.name
            elif r2.name == r.name: 
                consensus_seq = r2.sequence

        ### read coverage info from the cov gff for this sequence
        coverage_map={}
        for record in GffReader(self.args.proj+"/reseq/coverage_1.gff"):
            if record.seqid == name:
                for i in xrange(int(record.start), int(record.end)+1):
                    coverage_map[i] = record.getAttrVal('cov2') 

        ### read Qv info from the fastq for this sequence
        quality_map = {}
        i = 1
        for qual in r.quality:
            quality_map[i] = qual
            i += 1

        ### create annotation gff3 for this sequence
        ### comparing the same ref between the new and old profile allows us to keep track of the "canonical" coordinates
        gff3_annotation = create_annotation(comparison_seq, new_comparison_seq, consensus_seq, coverage_map, MSA_info, r.name, quality_map)     
        last_feature = None; feature_regions = {}
        for record in gff3_annotation:
            gDNA_annot_writer.writeRecord(record)       
            try:
                feature_regions[record.getAttrVal('feature')].append(record.getAttrVal('sequence_position'))
            except KeyError:
                feature_regions[record.getAttrVal('feature')] = [int(record.getAttrVal('sequence_position'))]
        for item in feature_regions.iteritems():
            feature = item[0]
            if 'exon' not in feature: continue
            bases = sorted(item[1])     
            record = Gff3Record(r.name) 
            record.type='region'
            record.start=bases[0]
            record.end=bases[-1]
            record.put('feature', feature )
            feature_writer.writeRecord(record)
            
        ### find best match among all profiles for our consensus sequence       
        ### write it out to a file
        best_gDNA_match, best_gDNA_match_score, gDNA_var_map = MSA_aligner(consensus_seq, output, r.name)           
        with open(self.args.proj+"/annotate/gDNA_allele_calls.txt", "a") as of:
            print>>of, "%s %s %s" % (r.name, best_gDNA_match, best_gDNA_match_score)    
        gDNA_calls[r.name] = locus
        
        ### create variant gff3 to document differences between consensus sequence and best matching reference
        gff3_var_annotation = create_var_annotation(consensus_seq, gDNA_var_map, r.name)    
        for record in gff3_var_annotation:
            gDNA_var_writer.writeRecord(record)

        ### create artificial cDNA sequence from the gDNA sequence
        ### create a dict between gDNA and cDNA position so that we can read in quality and coverage info to the cDNA gff
        ### and so we can keep track of where our artificial cDNA sequence fits in to the real sequence
        cDNA_consensus_sequence=''
        cDNA_base_counter=1
        gDNA_pos_to_cDNA_pos={}
        for i in xrange(len(gff3_annotation)):
            if gff3_annotation[i].type == 'base':       
                if gff3_annotation[i].getAttrVal('isexon'):
                    cDNA_consensus_sequence+=gff3_annotation[i].getAttrVal('basecall')
                    gDNA_pos_to_cDNA_pos[gff3_annotation[i].getAttrVal('sequence_position')] = cDNA_base_counter
                    cDNA_base_counter+=1
        cDNA_coverage_map={}
        cDNA_pos_to_gDNA_pos={}
        for item in gDNA_pos_to_cDNA_pos.iteritems():
            cDNA_pos_to_gDNA_pos[int(item[1])] = item[0]
        
        for item in coverage_map.iteritems():
            try:
                cDNA_coverage_map[gDNA_pos_to_cDNA_pos[item[0]]] = item[1]
            except:
                pass

        cDNA_quality_map={}
        for item in quality_map.iteritems():
            try:
                cDNA_quality_map[gDNA_pos_to_cDNA_pos[item[0]]] = item[1]
            except:
                pass

        with open(self.args.proj+"/annotate/resequenced_hap_con_cDNA.fasta", "a") as of:
            print>>of, ">"+r.name
            print>>of, cDNA_consensus_sequence

        ### now that we have established the cDNA sequence ... we start all over again with a cDNA MSA!
        ### get cDNA profile fn
        MSA_cDNA_info_fn = MSA_cDNA_info_fn_dict[locus]

        ### read in cDNA profile features
        MSA_cDNA_info = {}
        with open(MSA_cDNA_info_fn, "r") as f2:
            for line in f2:
                line = line.strip().split()
                MSA_cDNA_info[int(line[0])] = info._make( ( line[1], line[2], line[3] ) )
        MSA_cDNA_fn = MSA_cDNA_fn_dict[locus]
        with open(tmp_fn, "w") as of:
            print>>of, ">"+r.name       
            print>>of, cDNA_consensus_sequence

        ### add artificial cDNA sequence to the cDNA profile
        output = self.args.proj+"/annotate/"+name+"_cDNA.afa"
        if not os.path.isfile(output):
            try:
                muscle_output = check_output(". /mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh; \
                    muscle -profile -in1 %s -in2 %s -out %s 2> /dev/null" % (MSA_cDNA_fn, tmp_fn, output),
                    executable='/bin/bash', shell=True)
            except:
                log.info("MSA failed for ( %s )" % (name+"_cDNA") )
                os.remove(tmp_fn)
                continue

        ### read out a reference sequence from the profile
        f2 = FastaReader(MSA_cDNA_fn)
        letters_max = 0
        for r2 in f2:
            num_letters = ( r2.sequence.count('A') + r2.sequence.count('G') + r2.sequence.count('C') + r2.sequence.count('T') )
            if num_letters > letters_max:
                letters_max = num_letters
        f2 = FastaReader(MSA_cDNA_fn)
        for r2 in f2:
            num_letters = ( r2.sequence.count('A') + r2.sequence.count('G') + r2.sequence.count('C') + r2.sequence.count('T') )
            if num_letters == letters_max:
                comparison_seq = r2.sequence
                comparison_seq_name = r2.name
                break

        ### read out the same reference sequence from the NEW profile
        ### as well as our artificial cDNA consensus sequence
        f2 = FastaReader(output)
        for r2 in f2:
            if r2.name == comparison_seq_name:
                new_comparison_seq = r2.sequence
                new_comparison_seq_name = r2.name
            elif r2.name == r.name: 
                consensus_seq = r2.sequence

        ### create annotation gff3 for this sequence
        ### comparing the same ref between the new and old profile allows us to keep track of the "canonical" coordinates
        gff3_annotation = create_annotation(comparison_seq, new_comparison_seq, consensus_seq, cDNA_coverage_map, MSA_cDNA_info, r.name, cDNA_quality_map, coord_dict = cDNA_pos_to_gDNA_pos )  
        for record in gff3_annotation:
            cDNA_annot_writer.writeRecord(record)
        ## find best match among all profiles for our consensus sequence    
        best_cDNA_match, best_cDNA_match_score, cDNA_var_map = MSA_aligner(consensus_seq, output, r.name)           
        with open(self.args.proj+"/annotate/cDNA_allele_calls.txt", "a") as of:
            print>>of, "%s %s %s" % (r.name, best_cDNA_match, best_cDNA_match_score)
        cDNA_calls[r.name] = locus

        ### create variant gff3 to document differences between consensus sequence and best matching reference
        gff3_var_annotation = create_var_annotation(consensus_seq, cDNA_var_map, r.name, coord_dict = cDNA_pos_to_gDNA_pos )    
        for record in gff3_var_annotation:
            cDNA_var_writer.writeRecord(record)

        os.remove(tmp_fn)
    ### now we create a special gff that will allow us to do phase highlighting
    ### TODO: need to have a locus dict
    locus_counts = Counter(self.locus_dict.values())
    loci_to_compare=[]
    for item in locus_counts.iteritems():
        if item[1] == 2:
            loci_to_compare.append(item[0])

    phase_writer = GffWriter(self.args.proj+"/annotate/phase.gff")
    phase_writer.writeMetaData('pacbio-variant-version', '1.4')

    for locus in loci_to_compare:
        seqs_to_compare=[]
        for item in self.locus_dict.iteritems():
            if item[1] == locus:
                seqs_to_compare.append(item[0]+"|quiver")
        sequences = extract_sequence(self.args.proj+"/reseq/reseq_output_1.fasta", seqs_to_compare)
        write_fasta( sequences, tmp_fn, "w")
        ### we will align both seqs from the same locus to the MSA profile so that we can compare them
        output = self.args.proj+"/annotate/locus_"+locus+"_comparison.afa"
        MSA_fn = MSA_fn_dict[locus]
        if not os.path.isfile(output):
            try:
                muscle_output = check_output(". /mnt/secondary/Smrtanalysis/opt/smrtanalysis/etc/setup.sh; \
                    muscle -in %s -out %s 2> /dev/null" % (tmp_fn, output),
                    executable='/bin/bash', shell=True)
            except:
                log.info("MSA failed for ( %s )" % (locus) )
                os.remove(tmp_fn)
                continue
        sequences = extract_sequence(output, seqs_to_compare)
        if len(sequences) != 2: continue
        alignment = MSA_aligner(sequences[0].sequence, sequences[1].sequence, mode = 'string_to_string')
        qname = sequences[0].name
        tname = sequences[1].name
        q_vars = create_var_annotation(sequences[0].sequence, alignment.qvars, qname )          
        for record in q_vars:
            phase_writer.writeRecord(record)
        t_vars = create_var_annotation(sequences[1].sequence, alignment.tvars, tname )
        for record in t_vars:
            phase_writer.writeRecord(record)
    return

if __name__ == '__main__':
    import argparse

    desc = "A tool for annotating genomic HLA sequences"
    parser = argparse.ArgumentParser( description=desc )

    add = parser.add_argument
    add('stuff')

    args = parser.parse_args()

    log.basicConfig( level=logging.INFO )

