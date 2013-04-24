#!/home/UNIXHOME/jquinn/HGAP_env/bin/python

#################################################################################$$
# Copyright (c) 2011,2012, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without 
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright notice, this 
#   list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice, 
#   this list of conditions and the following disclaimer in the documentation 
#   and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its contributors 
#   may be used to endorse or promote products derived from this software 
#   without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED 
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS 
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#################################################################################$$


import os
import random
import sys

from BasH5Reader import BasH5Reader

MIN_LENGTH = 500
MIN_SCORE = 0.75

class BasH5Extractor( object ):
    """
    A tool for extracting sub-reads from BasH5 files
    """

    def __init__( self, filename,
                        min_length=None,
                        min_score=None,
                        dilution=None ):
        # First we convert the input file into a list of file-paths 
        self.file_list = self.parse_bash5_files( filename )
        # Next set remaining args or their defaults
        if min_length is None:
            self.min_read_length = MIN_LENGTH
        if min_score is None:
            self.min_read_score = MIN_SCORE
        self.dilution = dilution
        # Finally we validate all of the settings
        self.validate_settings()

    def validate_settings( self ):
        try:
            assert isinstance(self.min_read_length, int)
        except AssertionError:
            msg = "Minimum Length must be integer"
            raise TypeError( msg )
        try:
            assert isinstance(self.min_read_score, float)
        except AssertionError:
            msg = "Minimum Length must be integer"
            raise TypeError( msg )

    def parse_bash5_files( self, filename ):
        """Convert the input file to a list of absolute paths to Bash5s"""
        assert os.path.isfile( filename )
        if filename.endswith('.bas.h5'):
            return [os.path.abspath( filename )]
        elif filename.endswith('.fofn'):
            file_list = []
            with open( filename, 'r' ) as handle:
                for line in handle:
                    file_list.append( os.path.abspath(line.strip()) )
            return file_list
        else:
            msg = "Input file must be Bas.H5 of FOFN!"
            raise TypeError( msg )

    def extract_subreads( self ):
        for filename in self.file_list:
            for zmw in BasH5Reader( filename ):
                if zmw.readScore < self.min_read_score:
                    continue
                adapter_starts = [z.readStart for z in zmw.adapters]
                adapter_ends = [z.readEnd for z in zmw.adapters]
                for subread in zmw.subreads:
                    if len(subread.basecalls()) < self.min_read_length:
                        continue
                    if subread.readStart in adapter_ends and subread.readEnd in adapter_starts:
                        print ">fp_%d_%d_%d" % (subread.holeNumber, subread.readStart, subread.readEnd)
                        print subread.basecalls()
                    else:
                        print ">%d_%d_%d" % (subread.holeNumber, subread.readStart, subread.readEnd)
                        print subread.basecalls()
                        
    def __call__( self ):
        self.extract_subreads()

if __name__ == '__main__':
    extractor = BasH5Extractor( sys.argv[1] )
    extractor()
