import logging
import subprocess

PULSE_METRICS = "DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag"
NOISE_DATA = "-77.27,0.08654,0.00121"
MIN_ACCURACY = 0.75
MIN_LENGTH = 500
BESTN = 1
NCANDIDATES = 30

class SmrtAnalysisRunner( object ):
    "A tool for simplifying calls out to external SMRT Analysis tools"

    def __init__( self, setup, nproc=4 ):
        self.setup = setup
        self.nproc = nproc
        self.log = logging.getLogger(__name__)
        self.proc_num = 0

    def run_subprocess( self, proc_name, proc_args, proc_log=None, report=True ):
        # Log the process start and designated log file
        self.proc_num += 1
        proc_id = str(self.process).zfill(4)
        if proc_log is None:
            proc_log = 'smrt_analysis_%s.log' % proc_id
        self.log.info("Starting SMRT Analysis Process #%s (%s)..." % (proc_id, 
                                                                      proc_name))
        arg_string = ' '.join(proc_args)
        self.log.info('Running the following command:\n\t%s' % arg_string)
        self.log.info('Logging process to "%s"' % proc_log)
        # Run the specified process with the supplied args
        p = subprocess.Popen( proc_args, 
                              stderr=subprocess.STDOUT, 
                              stdout=proc_log )
        p.wait()
        # Read the process log into the global log
        if report:
            log_contents = self.read_log_file( proc_log )
            self.log.info("Process concluded with the following messages:\n%s\n" % log_contents)
        else:
            self.log.info("Process concluded\n")

    def read_log_file( self, log_file ):
        with open( log_file, 'r' ) as handle:
            lines = handle.readlines()
        return '\t'.join( lines )

    def check_output( self, output ):
        try:
            assert os.path.exists( output )
        except:
            msg = 'Expected output "%s" not found!' % output
            self.log.error( msg )
            raise IOError( msg )
        self.log.info('Expected output "%s" found' % output)

    def compare_sequences( self, query, reference, output):
        process_args = ['source', self.setup, ';',
                        'compareSequences.py',
                        '--info',
                        '--algorithm=blasr',
                        '-x', '-bestn', BESTN,
                        '-x', '-nCandidates', NCANDIDATES,
                        '--nproc=%s' % self.nproc,
                        '--placeRepeatsRandomly',
                        '--useQuality',
                        '--minAccuracy=%s' % MIN_ACCURACY,
                        '--minLength=%s' % MIN_LENGTH,
                        '--refSeqName=HLA',
                        '--noiseData=%s' % NOISE_DATA,
                        '--h5pbi', 
                        '--h5fn=%s' % output,
                        query, 
                        reference]
        self.run_subprocess( 'compareSequences.py', process_args )
        self.check_output( output )

    def load_pulses( self, basH5, cmpH5 ):
        process_args = ['source', self.setup, ';',
                        'loadPulses',
                        basH5,
                        cmpH5,
                        '-metrics', PULSE_METRICS]
        self.run_subprocess( 'loadPulses', process_args )
        self.check_output( cmpH5 )

    def sort_cmpH5( self, cmpH5 ):
        process_args = ['source', self.setup, ';',
                        'cmpH5tools.py', 'sort', cmpH5]
        self.run_subprocess( 'cmph5tools-sort', process_args )
        self.check_output( cmpH5 )

    def quiver( self, cmpH5, reference, fasta_out, fastq_out, gff_out ):
        process_args = ['source', self.setup, ';',
                        'variantCaller.py',
                        '--noEvidenceConsensusCall', 'reference',
                        '--algorithm=quiver',
                        '-j', self.nproc,
                        '-r', reference,
                        '-o', gff_out,
                        '-o', fasta_out,
                        '-o', fastq_out,
                        cmpH5]
        self.run_subprocess( 'Quiver', process_args, fasta_out )
        self.check_output( fasta_out )
        self.check_output( fastq_out )
        self.check_output( gff_out )

    def reference_uploader( self, fasta, name, output_dir ):
        process_args = ['source', self.setup, ';',
                        'referenceUploader', 
                        '-c', # Create new reference
                        '-f', fasta,
                        '-n', name,
                        '-p', output_dir]
        self.run_subprocess( 'referenceUploader', process_args )
        self.check_output( output_dir )

    def summarize_coverage( self, cmpH5, reference_dir, output ):
        process_args = ['source', self.setup, ';',
                        'summarizeCoverage.py',  cmpH5,
                        '--reference', reference_dir,
                        '--regionSize', REGION_SIZE] 
        self.run_subprocess( 'summarizeCoverage.py', process_args, 
                                                     output, 
                                                     report=False)
        self.check_output( output )
