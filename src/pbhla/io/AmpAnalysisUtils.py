from string import maketrans

from pbcore.io.base import ReaderBase, WriterBase, getFileHandle

from pbhla.io.AmpAnalysisIO import AmpAnalysis, AmpAnalysisReader, AmpAnalysisWriter
from pbhla.io.WhiteListReader import WhiteListReader
from pbhla.io.OrientationReader import OrientationReader

complement = maketrans('AGCTN', 'TCGAN')

def reverse_complement( record, length ):
    """
    Reverse-complement an AmpliconAnalysis record
    """
    return AmpAnalysis(BarcodeName   = record.BarcodeName,
                       FastaName     = record.FastaName,
                       CoarseCluster = record.CoarseCluster,
                       Phase         = record.Phase,
                       TotalCoverage = record.TotalCoverage,
                       Position      = str(length-int(record.Position)),
                       Base          = record.Base.translate( complement ),
                       QV            = record.QV,
                       Coverage      = record.Coverage,
                       HetScore      = record.HetScore)

def orient_amp_analysis( input_file, alignment_file, output_file=None ):
    """
    Reverse-complement, re-sort and output AmpliconAnalysis records
    """
    output_file = output_file or _get_output_file( input_file, 'oriented' )
    lengths = read_amp_analysis_lengths( input_file )
    orientations = OrientationReader( alignment_file )
    records = rev_comp_amp_analysis_records( input_file, orientations, lengths )
    # Re-sort and output the AA records now that they're in the correct orientation
    records = sorted(records, key=lambda x: int(x.Position))
    records = sorted(records, key=lambda x: x.FastaName)
    with AmpAnalysisWriter( output_file ) as writer:
        writer.write_header()
        for record in records:
            writer.write_record( record )
    return output_file

def subset_amp_analysis( input_file, whitelist_file, output_file=None ):
    """
    Subset the AmpliconAnalysis data from list of selected consensus sequences
    """
    output_file = output_file or _get_output_file( input_file, 'selected' )
    white_list = WhiteListReader( whitelist_file )
    with AmpAnalysisWriter( output_file ) as writer:
        writer.write_header()
        for record in AmpAnalysisReader( input_file ):
            if record.FastaName in white_list:
                writer.write_record( record )
    return output_file

def read_amp_analysis_lengths( input_file ):
    """
    Parse the lengths of each Amplicon Analysis consensus from their max positions
    """
    lengths = {}
    for record in AmpAnalysisReader( input_file ):
        try:
            lengths[record.FastaName] = max(int(record.Position), lengths[record.FastaName])
        except:
            lengths[record.FastaName] = int(record.Position)
    return lengths

def rev_comp_amp_analysis_records( input_file, orientations, lengths ):
    """
    Read the records from an AmpliconAnalysis file, reverse-complementing as needed
    """
    records = []
    for record in AmpAnalysisReader( input_file ):
        try:
            if orientations.is_reverse( record.FastaName ):
                records.append( reverse_complement( record, lengths[record.FastaName] ))
            else:
                records.append( record )
        except:
            records.append( record )
    return records

def _get_output_file( input_file, file_type ):
    basename = '.'.join( input_file.split('.')[:-1] )
    return '%s.%s.csv' % (basename, file_type)
