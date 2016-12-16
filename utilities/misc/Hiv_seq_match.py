#!/usr/bin/env python3

import sys
import re
import gzip
from optparse import OptionParser

def parse_options():
    parser = OptionParser( usage="%prog [options]", version="v1.0" )
    parser.add_option( '-f', '--file', dest='filename', help='fasta file with illumina reads containing viral fragements and LTR' )
    parser.add_option( '-o', '--output_token', dest='output_token', help="output token to output file with out an extension" )
    return parser

def main( filename, output_token ):
    #LTR = re.compile( 'GGAAAATCTCTAGCA' ) #HIV data LTR
    LTR = re.compile( 'CCGTAGTACTTCGGTACAACA' )  #REV sample LTR
    LINKER = re.compile( 'TAGTCCCTTAA' )   #Linker seq for both HIV and REV
    
    total_seq_count = 0
    seqs_with_LTR = 0
    seqs_with_LTR_long_reads = 0
    seqs_with_LTR_short_reads = 0
    seqs_with_LINKER = 0
    seqs_with_LINKER_long_reads = 0
    seqs_with_LINKER_short_reads = 0
    seqs_without_LINKER = 0
    
    LTR_and_LINKER_file = gzip.open( output_token + '.LTR_and_LINKER.fasta.gz', 'wb' )
    LTR_but_no_LINKER_file = gzip.open( output_token + '.LTR_but_no_LINKER.fasta.gz', 'wb' )
    LTR_but_no_LINKER_short_reads_file = gzip.open( output_token + '.LTR_but_no_LINKER_short_reads.fasta.gz', 'wb' )
    LTR_and_LINKER_short_reads_file = gzip.open( output_token + '.LTR_and_LINKER_short_reads.fasta.gz', 'wb' )
    
    with open( filename, 'r' ) as f:
        for line in f:
    
            header = line.strip()
            seq = f.readline().strip()

            total_seq_count += 1
            LTR_match = LTR.search( seq )

            if LTR_match:
                seqs_with_LTR += 1
                LINKER_match = LINKER.search( seq )
                if LINKER_match:
                    seqs_with_LINKER += 1
                    full_string = seq[ LTR_match.end() : LINKER_match.start() ]
                    if len( full_string ) >= 25 :
                        seqs_with_LTR_long_reads += 1
                        seqs_with_LINKER_long_reads += 1
                        LTR_and_LINKER_file.write( ( header + "\n" ).encode( 'utf-8' ) )
                        LTR_and_LINKER_file.write( ( full_string + "\n" ).encode( 'utf-8' ) )
                    else:
                        seqs_with_LTR_short_reads += 1
                        seqs_with_LINKER_short_reads += 1
                        LTR_and_LINKER_short_reads_file.write( ( header + "\n" ).encode( 'utf-8' ) )
                        LTR_and_LINKER_short_reads_file.write( ( full_string + "\n" ).encode( 'utf-8' ) )
                else:
                    seqs_without_LINKER += 1
                    full_string = seq[ LTR_match.end() : LTR_match.end() + 150 ]
                    if len( full_string ) >= 25 :
                        seqs_with_LTR_long_reads += 1
                        LTR_but_no_LINKER_file.write( ( header + "\n" ).encode( 'utf-8' ) )
                        LTR_but_no_LINKER_file.write( ( full_string + "\n" ).encode( 'utf-8' ) )
                    else:
                        seqs_with_LTR_short_reads += 1
                        LTR_but_no_LINKER_short_reads_file.write( ( header + "\n" ).encode( 'utf-8' ) )
                        LTR_but_no_LINKER_short_reads_file.write( ( full_string + "\n" ).encode( 'utf-8' ) )

            else:
                continue

    LTR_and_LINKER_file.close()
    LTR_but_no_LINKER_file.close()
    LTR_but_no_LINKER_short_reads_file.close()
    LTR_and_LINKER_short_reads_file.close()

    sys.stderr.write( "Total number of sequences: {:,}\n".format( total_seq_count ) )
    sys.stderr.write( "% of sequences with LTR: {:.2%} ({:,})\n".format( ( seqs_with_LTR / total_seq_count ), seqs_with_LTR ) )
    sys.stderr.write( "% of sequences with LTR (>=25bp): {:.2%} ({:,})\n".format( ( seqs_with_LTR_long_reads / total_seq_count ), seqs_with_LTR_long_reads ) )
    sys.stderr.write( "% of sequences with LTR (<25bp): {:.2%} ({:,})\n".format( ( seqs_with_LTR_short_reads /total_seq_count ), seqs_with_LTR_short_reads ) )
    sys.stderr.write( "% of sequences with LTR but with out LINKER: {:.2%} ({:,})\n".format( ( seqs_without_LINKER / total_seq_count ), seqs_without_LINKER ) )
    sys.stderr.write( "% of sequences with LTR and LINKER: {:.2%} ({:,})\n".format( ( seqs_with_LINKER / total_seq_count ), seqs_with_LINKER ) )
    sys.stderr.write( "% of sequences with LTR and LINKER (>=25bp): {:.2%} ({:,})\n".format( ( seqs_with_LINKER_long_reads / total_seq_count ), seqs_with_LINKER_long_reads ) )
    sys.stderr.write( "% of sequences with LTR and LINKER (<25bp): {:.2%} ({:,})\n".format( ( seqs_with_LINKER_short_reads / total_seq_count ), seqs_with_LINKER_short_reads ) )

if __name__ == '__main__':
    parser = parse_options( )
    ( options, args ) = parser.parse_args( )
    if options.filename == None:
        parser.error( "--file is a required argument" )
    if options.output_token == None:
        parser.error( "--output_token is a required argument" )
    main( options.filename, options.output_token )

sys.exit( 0 )
