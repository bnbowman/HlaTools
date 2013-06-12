def phase_subreads( self ):
    self.log.info("Starting the sub-read phasing process")
    # First 
    subread_files = []
    with open( self.subread_files ) as handle:
        for line in handle:
            info_tuple = file_info._make( line.strip().split() )
            subread_files.append( info_tuple )
    ### run through the subreads files to figure out which loci should actually be run through phasr
    ### if there are already two subread files for a single locus, theres little point in running through phasr
    ### as the two alleles are already well separated
    loci_count = Counter([ x.locus for x in subread_files ])
    to_be_phased={}
    for item in loci_count.iteritems():
        if item[1] == 1:
            to_be_phased[item[0]] = False
        else:
            to_be_phased[item[0]] = True
    
    seq_to_locus={}
    backbone_fasta = self.args.proj+"/phased/backbone.fa"
    self.hap_con_fasta = self.args.proj + "/phased/haplotype_consensus_sequences.fasta"
    for item in subread_files:
        ref_name = item.ref_name
        locus = item.locus
        fasta_file = item.fasta_fn
        ### if these are the unmapped reads we are going to run them through HGAp
        if locus == "unmapped":
            if self.args.HGAP:
                pass
                #phase_unmapped_reads()
            continue
        
        phasr_output = self.args.proj+"/phased/"+ref_name+"_haplotype_consensus.fa"

        ### extract the reference sequence that was used to "cluster" this group of reads
        sequence = extract_sequence(self.reference_sequences, [ ref_name ] )
        write_fasta(sequence, backbone_fasta, "w")
        
        try:
            assert os.path.isfile( fasta_file )
        except:
            msg = "Could not open ( %s ). Skipping." % fasta_file
            self.log.info( msg )
            raise IOError( msg )

        if os.path.isfile(phasr_output):
            for record in FastaReader(phasr_output):
                seq_to_locus[record.name] = locus
            self.log.info("Found phasr output for ( %s ) in ( %s ). Skipping." % (fasta_file, phasr_output))
            continue
        
        ### if this locus is not to be phased, we can just set max recursion level to 0, and phasr will build regular consensus
        if self.args.avoid_phasr and not to_be_phased[locus]:
            self.log.info('Running Gcon instead of Phasr for "%s"' % fasta_file)
            argstring = re.sub('--max_recursion_level [0-9]',
                               '--max_recursion_level 0',
                               self.phasr_argstring)
            argstring += ' --max_recursion_level 0'
        else:
            self.log.info('Running  Phasr on "%s"' % fasta_file)
            argstring = self.phasr_argstring
        phasr_cmd = "phasr %s --ref %s\
            --output %s \
            --log \
            %s" % ( fasta_file, backbone_fasta, phasr_output, self.phasr_argstring)
        ### run phasr
        self.log.info("phasr invoked: %s" % phasr_cmd)
        check_output(". /home/UNIXHOME/jquinn/HGAP_env/bin/activate; %s" % (phasr_cmd),  
                     executable='/bin/bash', shell=True)
        if os.path.isfile( phasr_output ):
            output_count = fasta_size(phasr_output)
            self.log.info('Phasr output %s seqs to "%s"' % ( output_count, phasr_output) )
            ### append results to multifasta
            runbash( "cat %s >> %s" % (phasr_output, self.hap_con_fasta) )
            for record in FastaReader(phasr_output):
                seq_to_locus[record.name] = locus
            ### create record of files created
            with open(self.args.proj+"/phased/phasr_output_files.txt", "a") as of:
                print >>of, "%s %s %s" % (phasr_output, ref_name, locus)
    # Finally, we write a summary of our Phasr outputs to file
    with open(self.args.proj+"/phased/phasr_output_seqs.txt", "w") as of:
        for item in seq_to_locus.iteritems():
            print >>of, "%s %s" % (item[0], item[1])
    phasr_output_count = fasta_size( self.hap_con_fasta )
    self.log.info("Phasr created %s sequence(s) total" % phasr_output_count )
    self.log.info("Phasing complete.\n")
