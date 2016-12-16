#!/usr/bin/env python3

"""
	Input: Reads barcodes from STDIN stream.
	Output: Returns 1 if there is no duplicate combination in barcodes
		Returns 0 if there is duplication
	Author: Mahesh Vangala
	Date:	July-29-2014
	
	Dependencies:	Requires Python3
				itertools
"""

import itertools
import sys
from optparse import OptionParser

def parse_options():
	parser = OptionParser( usage="%prog [options]", version="v1.0" )
	parser.add_option( '-f', '--file', dest='filename', help='filename with barcode information' )
	return parser

def main( file ):
	barcodes = []
	f = open( file, 'r' )
	header = f.readline()
	for line in f.readlines():
		line = line.strip()
		barcode = ( line.split( "," ) )[ 4 ]
		barcodes.append( barcode.strip() )
	f.close()
	global_pool = {}
	for barcode in barcodes:
		current_pool = {}
		for index in range( 0, len( barcode ) ):
			for base in ('A','G','C','T'):
				cur_bar = barcode[0:index] + base + barcode[index+1:]
				cur_bar = ''.join( cur_bar )
				current_pool[ cur_bar ] = barcode

		for key in current_pool.keys():
			if key in global_pool:
				sys.stderr.write( "MSG: barcode {0} occurred in more than one sample {1} and {2}".format( key, current_pool[key], global_pool[key] ) )
				return 0
			else:
				global_pool[ key ] = current_pool[ key ]
	sys.stderr.write( "MSG: No duplication found in samples" )
	return 1
				
if __name__ == '__main__':
	parser = parse_options( )
	(options, args) = parser.parse_args( )
	if options.filename == None :
		parser.error( "--file is a required argument" )
	dup_found = main( options.filename ) #returns 0 if duplicates found otherwise 1
	sys.stdout.write( str( dup_found ) )

sys.exit( 0 ) 
