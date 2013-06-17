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
import sys
import logging
import subprocess

from pbcore.io.FastqIO import FastqReader, FastqWriter
from pbcore.io.FastaIO import FastaRecord, FastaWriter  

PULSE_METRICS = 'DeletionQV,IPD,InsertionQV,PulseWidth,QualityValue,MergeQV,SubstitutionQV,DeletionTag'

log = logging.getLogger()

class ClusterResequencer(object):
    """
    A tool for resequencing clusters of rDNA sequences    
    """

    def __init__(self, reads, reference, fofn, setup=None, output='Resequencing', nproc=1):
        self.reads = reads
        self.reference = reference
        self.fofn = fofn
        self.setup = setup
        self.output = output
        self.nproc = nproc
        self.parseSequenceData()

    def validateSettings(self):
        self.filterPlsH5 = which('filterPlsH5.py')
        self.compareSequences = which('compareSequences.py')
        self.cmph5tools = which('cmph5tools.py')
        self.loadPulses = which('loadPulses')
        self.variantCaller = which('variantCaller.py')

    ####################
    # Instance Methods #
    ####################

    def parseDistances(self):
        distances = []
        with open( self.listFile, 'r' ) as handle:
            for line in handle:
                parts = line.split()
                distance = self.convertDistance( parts[0] )
                distances.append( distance )
        return distances

    def parseSequenceData(self):
        self.sequenceData = {}
        for fastqRecord in FastqReader( self.ccsFile ):
            zmw = getZmw( fastqRecord.name )
            self.sequenceData[zmw] = fastqRecord

    def trimClusterNames(self, clusters):
        trimmed = []
        for cluster in clusters:
            cluster = [getZmw(c) for c in cluster]
            trimmed.append( frozenset(cluster) )
        return trimmed

    def getClusterReads(self, cluster):
        reads = []
        for ccsZmw in cluster:
            try:
                ccsRead = self.sequenceData[ccsZmw]
            except KeyError:
                #raise Warning("No CCS read found for '%s', skipping..." % ccsZmw)
                continue
            reads.append( ccsRead )
        return reads

    def findLongestRead(self, reads):
        lengths = [len(read.sequence) for read in reads]
        maxLength = max(lengths)
        longestReads = [read for read in reads
                             if len(read.sequence) == maxLength]
        return longestReads[0]

    def outputClusterWhitelist(self, cluster, count):
        print "Creating Whitelist for Cluster #%s" % count
        whiteListFile = 'cluster%s_whitelist.txt' % count
        if os.path.exists( whiteListFile ):
            return whiteListFile
        with open( whiteListFile, 'w' ) as handle:
            for zmw in cluster:
                handle.write(zmw + '\n')
        return whiteListFile

    def outputClusterReference(self, reference, count):
        print "Creating reference sequence for Cluster #%s" % count
        referenceFile = 'cluster%s_reference.fasta' % count
        if os.path.exists( referenceFile ):
            return referenceFile
        # Rename the "Reference" sequence to the cluster
        referenceFasta = FastaRecord("Cluster%s" % count,
                                     reference.sequence)
        with FastaWriter( referenceFile ) as handle:
            handle.writeRecord( referenceFasta )
        return referenceFile

    def outputRepresentativeRead(self, representativeRead, count):
        print "Creating representative sequence file Cluster #%s" % count
        representativeFile = 'cluster%s_represent.fastq' % count
        if os.path.exists( representativeFile ):
            return representativeFile
        with FastqWriter( representativeFile ) as handle:
            handle.writeRecord( representativeRead )
        return representativeFile

    def createRgnH5(self, whiteListFile, count):
        print "Creating Region Table for Cluster #%s" % count
        outputDir = 'cluster%s_regionTables' % count
        outputFofn = 'cluster%s_regionTables.fofn' % count
        if os.path.exists( outputDir ) and os.path.exists( outputFofn ):
            return outputFofn
        outputDirArg = '--outputDir=%s' % outputDir
        outputFofnArg = '--outputFofn=%s' % outputFofn
        filterArg = '--filter=ReadWhitelist=%s,MinReadScore=0.75' % whiteListFile
        p = subprocess.Popen( [self.filterPlsH5, 
                               self.sequenceFile,
                               outputDirArg, 
                               outputFofnArg,
                               filterArg] )
        p.wait()
        print "Region Table Created Successfully"
        return outputFofn

    def createCmpH5(self, referenceFile, rgnH5File, count):
        print "Creating a CMP.H5 for Cluster #%s" % count
        cmpH5File = 'cluster%s.cmp.h5' % count
        if os.path.exists( cmpH5File ):
            return cmpH5File
        p = subprocess.Popen( [self.compareSequences, 
                               '--minAccuracy=0.75',
                               '--minLength=500',
                               '--useQuality',
                               '--h5pbi',
                               '--info',
                               '--nproc=4',
                               '-x', '-bestn', '1',
                               '--nproc=%s' % self.numProc,
                               '--regionTable=%s' % rgnH5File,
                               '--algorithm=blasr',
                               '--noiseData=-77.27,0.08654,0.00121',
                               '--h5fn=%s' % cmpH5File,
                               self.sequenceFile,
                               referenceFile] )
        p.wait()
        return cmpH5File

    def sortCmpH5(self, cmph5File, count):
        print "Sorting the CmpH5 for Cluster #%s" % count
        sortedCmpH5File = 'cluster%s.sorted.cmp.h5' % count
        if os.path.exists( sortedCmpH5File ):
            return sortedCmpH5File
        p = subprocess.Popen( [self.cmph5tools,
                               'sort',
                               '--outFile=%s' % sortedCmpH5File,
                               cmph5File] )
        p.wait()
        return sortedCmpH5File

    def loadPulsesIntoCmpH5(self, sortedCmpH5File, count):
        print "Loading pulse data into the CmpH5 for Cluster #%s" % count
        p = subprocess.Popen( [self.loadPulses,
                               self.sequenceFile,
                               sortedCmpH5File,
                               '-metrics',
                               PULSE_METRICS] )
        return True

    def runQuiver(self, referenceFile, sortedCmpH5File, count):
        print "Running Quiver-consensus on Cluster #%s" % count
        consensusFile = 'cluster%s_consensus.fastq' % count
        if os.path.exists( consensusFile ):
            return consensusFile
        p = subprocess.Popen( [self.variantCaller,
                               '--algorithm=quiver',
                               '--numWorkers=%s' % self.numProc,
                               '--reference=%s' % referenceFile,
                               '--outputFile=%s' % consensusFile,
                               sortedCmpH5File])
        p.wait()
        return consensusFile

    def combineOutputSequences(self, outputSequenceFiles):
        print "Combining Consensus and Representative sequences"
        outputSequences = []
        for filename in outputSequenceFiles:
            for record in FastqReader( filename ):
                outputSequences.append( record )
        return outputSequences

    def outputCombinedSequences(self, combinedSequences ):
        print "Writing finished sequences to file"
        with FastqWriter( self.outputFile ) as handle:
            for record in combinedSequences:
                handle.writeRecord( record )

    def __call__(self):
        outputSequenceFiles = []
        # Select the appropriate distance, and parse the matching clusters
        distances = self.parseDistances()
        trimmedClusters = self.trimClusterNames( clusters )
        # Iterate over the clusters, generating consensuses
        for count, cluster in enumerate( trimmedClusters ):
            count = str(count+1).zfill(4)
            print "Analyzing cluster #%s now..." % (count)
            reads = self.getClusterReads( cluster )
            # If we have enought reads
            if len(reads) >= self.minClusterSize:
                print "%s ZMWs found (of %s), generating consensus..." % \
                                                 (len(reads), len(cluster))
                # First we select the longest CCS read from the cluster
                longest = self.findLongestRead( reads )
                referenceFile = self.outputClusterReference( longest, count )
                # Second we create a Rng.H5 file to mask other reads from Blasr
                whiteListFile = self.outputClusterWhitelist( cluster, count )
                rgnH5File = self.createRgnH5( whiteListFile, count )
                # Third we create a sorted CmpH5
                cmpH5File = self.createCmpH5( referenceFile, rgnH5File, count )
                sortedCmpH5File = self.sortCmpH5( cmpH5File, count )
                # Fourth we load rich QV data and run Quiver
                self.loadPulsesIntoCmpH5( sortedCmpH5File, count )
                consensusFile = self.runQuiver( referenceFile, 
                                                sortedCmpH5File,
                                                count )
                # Finally we record the name of the output file
                outputSequenceFiles.append( consensusFile )
            # Otherwise, we select one "Best" read to represent the cluster
            else:
                print "%s ZMWs found (of %s), skipping consensus..." % \
                                               (len(reads), len(cluster))
                reads = self.getClusterReads( cluster )
                representRead = reads[0]
                representFile = self.outputRepresentativeRead( representRead, 
                                                               count )
                outputSequenceFiles.append( representFile )
        # Finally we combine and trim all of the output Files
        combinedSequences = self.combineOutputSequences( outputSequenceFiles )
        self.outputCombinedSequences( combinedSequences )

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    add = parser.add_argument
    add('read_file', metavar='READS', 
        help="Fasta or Fastq file of CCS sequences")
    add('ref_file',  metavar='REF',  
        help="Mothur list file of cluster data")
    add('fofn_file', metavar='FOFN', 
        help="BasH5 or FOFN of sequence data")
    add('--setup', metavar='SETUP',
        help='SMRT Analysis setup script')
    add('--output', default='reseq', metavar='DIR',
        help="Specify a directory for intermediate files")
    add('--nproc', type=int. default=1, metavar='INT',
        help="Number of processors to use")
    args = parser.parse_args()

    resequencer = rDnaResequencer(args.read_file, 
                                  args.ref_file, 
                                  args.fofn_file, 
                                  args.setup,
                                  args.output, 
                                  args.nproc)
    resequencer()
