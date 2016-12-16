#!/usr/bin/env python

import HTSeq
import numpy
import sys
import itertools

sam_file = HTSeq.SAM_Reader( sys.argv[1] )
genes = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
gtf_file = HTSeq.GFF_Reader( sys.argv[2], end_included=True )

counts = {}
for feature in gtf_file:
    if feature.type == "gene":
        genes[ feature.iv ] += feature.name
        counts[ feature.name ] = 0

for aln in sam_file:
    if aln.aligned:
        for iv2, step_set in genes[ aln.iv ].steps( ):
            for gene in list( step_set ):
                    counts[ gene ] += 1


for gene in sorted( counts.keys( ) ):
    if counts[ gene ] > 0:
        print "\t".join( ( gene, str( counts[ gene ] ) ) )
