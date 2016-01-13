#!/bin/bash

#input files
INPUTDIR=input/Bouchenak/
TREE=$INPUTDIR/bestTree.phy
ALN=$INPUTDIR/aln.exported.nex 

OUTDIR=output/
mkdir -p $OUTDIR

PERLSCRIPTDIR=perlScripts/
if [ -d "$PERLSCRIPTDIR" ];then
    $PERLSCRIPTDIR/parseNexusCharSetToDAV.pl -f $ALN > $OUTDIR/perlSubsets || exit
    $PERLSCRIPTDIR/displaysub.pl $TREE $OUTDIR/perlSubsets > $OUTDIR/perlSubtrees || exit
    $PERLSCRIPTDIR/maketriplets.pl $OUTDIR/perlSubtrees > $OUTDIR/perlTriplets || exit
    $PERLSCRIPTDIR/countParents.pl $OUTDIR/perlTriplets > $OUTDIR/perlParentCount || exit
    #$PERLSCRIPTDIR/sameTerrace.pl $INPUTDIR/100pars.phy $OUTDIR/perlSubsets > $OUTDIR/perlTerraceList || exit
fi

TERRAPHYDIR=../scripts/
MAIN=$TERRAPHYDIR/terraphy.main.py

#making the subsets currently takes a while in terraphy, so can copy over 
#the perl version if desired
#cp $OUTDIR/perlSubsets $OUTDIR/subsets

if [ ! -e $OUTDIR/subsets ];then
    echo SUBSETS
    $MAIN --coverage --alignment-file $ALN > $OUTDIR/subsets || exit
fi

rm -f $OUTDIR/subtrees $OUTDIR/triplets $OUTDIR/parents $OUTDIR/build.tre $OUTDIR/strict.tre

echo SUBTREES
$MAIN --display --parent-tree $TREE --subset-file $OUTDIR/subsets > $OUTDIR/subtrees || exit

echo TRIPLETS
$MAIN -t --subtree-file $OUTDIR/subtrees > $OUTDIR/triplets || exit
#grep -v "[ ].*[ ].*[ ].*" temp >> $OUTDIR/triplets || exit
#some giberish sometimes gets output to the stream
perl -p -i -e 's/^.*1034h//g' $OUTDIR/triplets

echo PARENTS
$MAIN --parents --triplet-file $OUTDIR/triplets > $OUTDIR/parentCount || exit

echo BUILD
$MAIN --build --triplet-file $OUTDIR/triplets > $OUTDIR/build.tre || exit

echo BUILD
$MAIN --strict --triplet-file $OUTDIR/triplets > $OUTDIR/strict.tre || exit

echo TERRACES
../scripts/terraphy.main.py --list-terraces --subset-file $OUTDIR/subsets --treefiles-to-assign $INPUTDIR/100pars.tre > $OUTDIR/terraceList

