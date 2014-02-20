#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import logging

from pbcore.io.BasH5IO import BasH5Collection, BasH5Reader, BaxH5Reader
from pbhla.utils import validate_file

log = logging.getLogger()

def get_bash5_reader( bash5_file ):
    bash5_file = validate_file( bash5_file )
    if bash5_file.endswith('.fofn'):
        c = BasH5Collection( bash5_file )
        assert len(c.movieNames) == 1       # Currently we only accept single-movies
        return c
    elif bash5_file.endswith('.bas.h5'):
        return BasH5Reader( bash5_file )
    elif bash5_file.endswith('.bax.h5'):
        return BaxH5Reader( bash5_file )
    else:
        raise ValueError

def get_zmw( bash5, zmw_num ):
    if isinstance( bash5, BasH5Reader ) or \
       isinstance( bash5, BaxH5Reader ):
        return bash5[zmw_num]
    elif isinstance( bash5, BasH5Collection ):
        movie = bash5.movieNames[0]
        zmw_name = '{0}/{1}'.format(movie, zmw_num)
        return bash5[zmw_name]
    else:
        raise ValueError

def filter_zmw_list( bash5, zmw_list, min_snr=None ):
    filtered_list = []
    for zmw_num in zmw_list:
        zmw = get_zmw( bash5, zmw_num )
        if min_snr and min( zmw.zmwMetric("HQRegionSNR")) < min_snr:
            continue
        filtered_list.append( zmw_num )
    fraction = 100 * round(len(filtered_list)/float(len(zmw_list)), 3)
    log.info("{0} of {1} ZMWs ({2}%) passed all quality filters".format(len(filtered_list),
                                                                         len(zmw_list),
                                                                         fraction))
    return filtered_list