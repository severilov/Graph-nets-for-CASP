#!/bin/bash

VORONOTADIR=../../voronota
MODELSDIR=../../data/CASP9/models
TARGETSDIR=../../data/CASP9/targets
OUTPUTDIR=../cad_results/CASP9

for INDIR in $MODELSDIR/*
do
	INDIRBASENAME=$(basename $INDIR .pdb)
	for INFILE in $INDIR/*
	do
		$VORONOTADIR/voronota-cadscore \
	  	-t $TARGETSDIR/${INDIRBASENAME}.pdb \
	  	-m $INFILE
	done > $OUTPUTDIR/${INDIRBASENAME}_cad_scores
done
