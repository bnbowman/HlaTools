import logging, subprocess, os

SIMILARITY = 0.95

log = logging.getLogger()

def cd_hit_est(query, output, similarity=SIMILARITY):
    similarity = str(similarity)

    log.info("Running Cd-Hit-Est with the following options:")
    log.info("\tQuery: " + query)
    log.info("\tOutput: " + output)
    log.info("\tSimilarity: " + similarity)

    log.info("Executing the Cd-Hit-Est command")
    command = ['cd-hit-est', '-i', query, 
                             '-o', output,
                             '-c', similarity]
    subprocess.check_call( command, stdout=None, stderr=None )
    log.info("Cd-Hit-Est filtering complete")

    # Remove unused files
    clean_up( output )

def clean_up( output ):
    log.info("Cleaning up unused Cd-Hit-Est output files")
    for ext in ['.clstr', '.bak.clstr']:
        filename = output + ext
        try:
            os.remove( filename )
        except:
            pass
