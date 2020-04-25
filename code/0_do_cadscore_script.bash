#!/bin/bash

VORONOTADIR=../../voronota
INPUTDIR=../../data/CASP12
OUTPUTDIR=../cad_results

for INDIR in $INPUTDIR/*
do
	INDIRBASENAME=$(basename $INDIR .pdb)
	for INFILE in $INDIR/*
	do
		$VORONOTADIR/voronota-cadscore \
	  	-t $INDIR/${INDIRBASENAME}.pdb \
	  	-m $INFILE
	done > $OUTPUTDIR/${INDIRBASENAME}_cad_scores
done
