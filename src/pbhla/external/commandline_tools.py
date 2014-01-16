import os, logging, subprocess

from pbhla.utils import check_output_file

log = logging.getLogger()

def run_blasr(query, reference, args, verbose=False):
    command_args = create_blasr_command(query, reference, args)
    if verbose:
        log_command( command_args )
    execute_command( command_args )
    if 'out' in command_args:
        check_output_file( command_args['out'] )

def run_muscle( args ):
    command_args = create_muscle_command( args )
    log_command( command_args )
    execute_command( command_args )

def run_hmmsearch(query, reference, args):
    command_args = create_hmmsearch_command(query, reference, args)
    log_command( command_args )
    execute_command( command_args )
    if 'domtblout' in command_args:
        check_output_file( command_args['domtblout'] )

def create_muscle_command( args ):
    log.debug("Converting supplied arguments into a Blasr commandline")
    command_args = ['muscle']
    for arg, value in args.iteritems():
        arg = '-' + str(arg)
        if value is True:
            command_args += [arg]
        else:
            command_args += [arg, str(value)]
    return command_args

def create_blasr_command(query, reference, args):
    log.debug("Converting supplied arguments into a Blasr commandline")
    command_args = ['blasr', query, reference]
    for arg, value in args.iteritems():
        arg = '-' + str(arg)
        if value is True:
            command_args += [arg]
        else:
            command_args += [arg, str(value)]
    return command_args

def create_hmmsearch_command(query, reference, args):
    log.debug("Converting supplied arguments into an Hmmsearch commandline")
    command_args = ['hmmsearch']
    for arg, value in args.iteritems():
        arg = '--' + str(arg)
        if value is True:
            command_args += [arg]
        else:
            command_args += [arg, str(value)]
    command_args += [query, reference]
    return command_args

def log_command( args ):
    args = list( args )
    command = args.pop(0).capitalize()
    log.debug('Executing "%s" with the following options:' % command)

    # Log any positional arguments
    while not args[0].startswith('-'):
        arg = args.pop(0)
        log.debug('\tArgument: %s' % arg)

    # Log any remaining options
    while True:
        if len(args) == 0:
            break
        option_parts = [args.pop(0)]
        while len(args) and not args[0].startswith('-'):
            part = args.pop(0)
            option_parts.append( part )
        log.debug('\tOption: "%s"' % ' '.join(option_parts))

def execute_command( command_args ):
    command = command_args[0].capitalize()
    log.debug('Executing "%s" command as subprocess' % command)
    with open('/dev/null', 'w') as handle:
        subprocess.check_call( command_args,
                               stdout=handle,
                               stderr=subprocess.STDOUT )
    log.debug("Subprocess finished successfully")
