import logging, subprocess

class BlasrRunner( object ):

    def __init__(self, query, reference, output, args={}):
        self.query = query
        self.reference = reference
        self.output = output
        self.blasr_args = args
        self.initialize_logging()
        self.run()

    def initialize_logging(self):
        self.log = logging.getLogger(__name__)
        self.log.info("Initializing BlasrRunner with the following:")
        self.log.info("\tQuery: {0}".format(self.query))
        self.log.info("\tReference: {0}".format(self.reference))
        self.log.info("\tOutput: {0}".format(self.output))
        for key, value in self.blasr_args.iteritems():
            self.log.info("\t{0}: {1}".format(key, value))

    def run(self):
        self.create_command()
        self.run_command()

    def create_command(self):
        self.log.info("Converting supplied arguments into a Blasr commandline")
        self.command_args = ['blasr', self.query, self.reference, '-out', self.output]
        for arg, value in self.blasr_args.iteritems():
            if arg == 'output_type':
                if value == 'sam':
                    self.command_args.append('-sam')
                elif value == 'm1':
                    self.command_args += ['-m', '1']
                elif value == 'm4':
                    self.command_args += ['-m', '4']
                elif value == 'm5':
                    self.command_args += ['-m', '5']
                else:
                    msg = 'Unrecognized output type! "{0}"'.format(value)
            else:
                self.command_args += ['-{0}'.format(arg), str(value)]

    def run_command(self):
        self.log.info("Executing the Blasr command")
        subprocess.check_call( self.command_args )
        self.log.info("Blasr alignment finished")
