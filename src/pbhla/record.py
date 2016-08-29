#! /usr/bin/env python

__author__ = 'Brett Bowman'
__email__ = 'bbowman@pacificbiosciences.com'

import re
import logging

log = logging.getLogger()

def consensus_size( record ):
    assert hasattr( record, 'name' )
    if re.search( 'NumReads\d+$', record.name.strip() ):
        return int( record.name.strip().split('NumReads')[1] )
    else:
        msg = 'Record has no "NumReads" attribute to return'
        log.error( msg )
        raise ValueError( msg )

def record_accuracy( record, precision=7 ):
    assert hasattr( record, 'quality' )
    p_values = [quality_to_p(qv) for qv in record.quality]
    average = sum(p_values)/len(p_values)
    return round(average, precision)

def quality_to_p(qv):
    return 1-(10**(-1*float(qv)/10.0))
