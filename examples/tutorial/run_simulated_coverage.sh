#!/bin/bash

#For the purposes of this script you can specify which terraphy program to use
#remove the leading # to uncomment one of the below TERRAPHY= as needed

#this requires that Terraphy is globally installed
TERRAPHY=terraphy.main.py
#if you have a 64 bit mac this shgould work, and won't even require a python installation
#TERRAPHY=../scripts/terraphy.main.standalone.OSX

#This script runs through the typical preprocessing steps 
#and many of the analyses implemented by the software

#INPUT files: parent tree and alignment
INDIR=input
TREE=$INDIR/Springer.catarrhiniNoMacaca.prunedOrig.tre

OUTDIR=simulatedCoverage.output
mkdir -p $OUTDIR

#PREVIOUS OUTPUT WILL BE OVERWRITTEN EACH TIME THIS SCRIPT IS CALLED
#SEE the examples/test.sh SCRIPT FOR AN EXAMPLE OF HOW TO NOT RECOMPUTE EXISTING FILES

#the intermediate files created by preproccesing of the input tree and alignment
#are named "subsets", "subtrees" and "triplets"
#and are used as input for the main functions of the software.  They will be written to the directory "preprocessing"
PREDIR=simulatedCoverage.preprocessing

for COV in  0.2 0.3 0.5 
do
NAME=coverage$COV
SUBSETS=$PREDIR/subsets.$NAME
SUBTREES=$PREDIR/subtrees.$NAME
TRIPLETS=$PREDIR/triplets.$NAME

#skip this, since using simulated subset files
#echo CREATING SUBSETS
#the subsets file requires the alignement as input 
#$TERRAPHY --coverage --alignment-file $ALN > $SUBSETS || exit

echo CREATING SUBTREES
#the subtrees file requires the parent tree and subsets file as input
$TERRAPHY --display --parent-tree-file $TREE --subset-file $SUBSETS > $SUBTREES || exit

echo CREATING TRIPLETS
#the triplets file requires the subtrees file  as input
$TERRAPHY -t --subtree-file $SUBTREES > $TRIPLETS || exit

PARENTS=$OUTDIR/parentCount.$NAME
echo COUNTING TREES ON TERRACE, i.e. "parent trees"
$TERRAPHY --count-parents --triplet-file $TRIPLETS > $PARENTS || exit

BUILD=$OUTDIR/build.$NAME.tre 
echo MAKING BUILD CONSENSUS TREE OF TERRACE TREES
$TERRAPHY --build  --annotate-clades --triplet-file $TRIPLETS > $BUILD || exit

PARENTS=$OUTDIR/allParents.$NAME.tre
#this function is generally impractical on most datasets and is commented out here
#echo GENERATING ALL TREES ON TERRACE
#$TERRAPHY --generate-parents --triplet-file $TRIPLETS > $PARENTS || exit

STRICT=$OUTDIR/strict.$NAME.tre 
#this function is impractical on many datasets and is commented out here
echo MAKING STRICT CONSENSUS OF TREES ON TERRACE
$TERRAPHY --strict --annotate-clades --triplet-file $TRIPLETS > $STRICT || exit

LIST=$OUTDIR/terraceList.$NAME
INPUTTREES=$INDIR/21LikelihoodTrees.tre
echo ASSIGNING TREES TO TERRACES
#This is a optional specialized functionality to determine which of a set of input trees lie on the same terrace
#$TERRAPHY --list-terraces --subset-file $SUBSETS --treefiles-to-assign $INPUTTREES > $LIST || exit

done
