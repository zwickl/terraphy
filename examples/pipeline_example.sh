#!/bin/bash

#This script runs through the typical preprocessing steps 
#and many of the analyses implemented by the software

#input files
#This example uses data from the Bouchenak dataset, found in the input/Bouchenak directory
#minimum input is a "parent" tree from the terrace to be analyzed, and a NEXUS format alignment 
#file containing a data matrix and a sets block with charsets that define the boundaries of the
#genes, loci or "partitions" analyzed in the partitioned analysis
DATA=Bouchenak
INPUTDIR=input/$DATA
TREE=$INPUTDIR/bestTree.phy
ALN=$INPUTDIR/aln.exported.nex 

#output files will be placed in the pipeline.output directory
OUTDIR=pipeline.output 
mkdir -p $OUTDIR

#the terraphy script in the in the current directory tree is used in this script
TERRAPHYDIR=../scripts/
MAIN=$TERRAPHYDIR/terraphy.main.py

#the intermediate files created by preproccesing of the input tree and alignment
#are named "subsets", "subtrees" and "triplets"
#which  are used as input for the main functions of the software
SUBSETS=$OUTDIR/subsets
SUBTREES=$OUTDIR/subtrees
TRIPLETS=$OUTDIR/triplets

echo CREATING SUBSETS
#the subsets file requires the alignement as input 
$MAIN --coverage --alignment-file $ALN > $SUBSETS || exit

echo CREATING SUBTREES
#the subtrees file requires the parent tree and subsets file as input
$MAIN --display --parent-tree-file $TREE --subset-file $SUBSETS > $SUBTREES || exit

echo CREATING TRIPLETS
#the triplets file requires the subtrees file  as input
$MAIN -t --subtree-file $SUBTREES > $TRIPLETS || exit
#grep -v "[ ].*[ ].*[ ].*" temp >> $TRIPLETS || exit
#some giberish sometimes gets output to the stream
#perl -p -i -e 's/^.*1034h//g' $TRIPLETS

PARENTS=$OUTDIR/parentCount 
echo COUNTING TREES ON TERRACE, i.e. "parent trees"
$MAIN --count-parents --triplet-file $TRIPLETS > $PARENTS || exit

BUILD=$OUTDIR/build.tre 
echo MAKING BUILD CONSENSUS TREE OF TERRACE TREES
$MAIN --build  --triplet-file $TRIPLETS > $BUILD || exit

PARENTS=$OUTDIR/allParents.tre
#this function is generally impractical on most datasets and is commented out here
#echo GENERATING ALL TREES ON TERRACE
#$MAIN --generate-parents --triplet-file $TRIPLETS > $PARENTS || exit

SAMP=100
SAMPLES=$OUTDIR/samp"$SAMP".tre 
echo SAMPLING $SAMP TREES FROM TERRACE
$MAIN --sample-parents $SAMP --triplet-file $TRIPLETS > $SAMPLES || exit

STRICT=$OUTDIR/strict.tre 
#this function is impractical on many datasets and is commented out here
#echo MAKING STRICT CONSENSUS OF TREES ON TERRACE
#$MAIN --strict --annotate-clades --triplet-file $TRIPLETS > $STRICT || exit

LIST=$OUTDIR/terraceList
INPUTTREES=$INPUTDIR/100pars.tre
echo ASSIGNING TREES TO TERRACES
#This is a optional specialized functionality to test whether a set of input trees lies on a single or
#multiple terraces
../scripts/terraphy.main.py --list-terraces --subset-file $SUBSETS --treefiles-to-assign $INPUTTREES > $LIST


