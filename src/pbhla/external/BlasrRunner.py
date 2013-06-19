import os, logging, subprocess

class BlasrRunner( object ):

    def __init__(self, query, reference, output, args={}):
        self.query = query
        self.reference = reference
        self.output = output
        self.output_type = output.split('.')[-1]
        self.blasr_args = args
        self.initialize_logging()
        self.run()

    def initialize_logging(self):
        self.log = logging.getLogger(__name__)
        self.log.info("Initializing BlasrRunner with the following:")
        self.log.info("\tQuery: {0}".format(os.path.basename(self.query)))
        self.log.info("\tReference: {0}".format(os.path.basename(self.reference)))
        self.log.info("\tOutput: {0}".format(self.output))
        for key, value in self.blasr_args.iteritems():
            self.log.info("\t{0}: {1}".format(key, value))

    def run(self):
        self.create_command()
        self.run_command()

    def create_command(self):
        self.log.info("Converting supplied arguments into a Blasr commandline")
        self.command_args = ['blasr', self.query, self.reference, 
                             '-out', self.output,
                             '-noSplitSubreads']
        for arg, value in self.blasr_args.iteritems():
            self.command_args += ['-{0}'.format(arg), str(value)]
        # Add a command for the output filetype 
        if self.output_type == 'sam':
            self.command_args += ['-sam']
        elif self.output_type == 'm1':
            self.command_args += ['-m', '1']
        elif self.output_type == 'm4':
            self.command_args += ['-m', '4']
        elif self.output_type == 'm5':
            self.command_args += ['-m', '5']
        else:
            msg = 'Unrecognized output type! "{0}"'.format(self.output_type)
            self.log.info( msg )
            raise TypeError( msg )

    def run_command(self):
        self.log.info("Executing the Blasr command")
        subprocess.check_call( self.command_args )
        self.log.info("Blasr alignment finished")
