#!/bin/bash

#DATA=Pyron
DATA=Bouchenak

#input files
INPUTDIR=input/$DATA
TREE=$INPUTDIR/bestTree.phy
ALN=$INPUTDIR/aln.exported.nex 

#OUTDIR=output.$DATA
OUTDIR=output.$DATA.test
mkdir -p $OUTDIR

PERLSCRIPTDIR=_perlScripts/

if [ -d "$PERLSCRIPTDIR" ];then
    $PERLSCRIPTDIR/parseNexusCharSetToDAV.pl -f $ALN > $OUTDIR/perlSubsets || exit
    $PERLSCRIPTDIR/displaysub.pl $TREE $OUTDIR/perlSubsets > $OUTDIR/perlSubtrees || exit
    $PERLSCRIPTDIR/maketriplets.pl $OUTDIR/perlSubtrees > $OUTDIR/perlTriplets || exit
    $PERLSCRIPTDIR/countParents.pl $OUTDIR/perlTriplets > $OUTDIR/perlParentCount || exit
    #$PERLSCRIPTDIR/sameTerrace.pl $INPUTDIR/100pars.phy $OUTDIR/perlSubsets > $OUTDIR/perlTerraceList || exit
fi

TERRAPHYDIR=../scripts/
MAIN=$TERRAPHYDIR/terraphy.main.py


SUBSETS=$OUTDIR/subsets
SUBTREES=$OUTDIR/subtrees
TRIPLETS=$OUTDIR/triplets

#making the subsets currently takes a while in terraphy, so can copy over 
#the perl version if desired
#cp $OUTDIR/perlSubsets $SUBSETS

if [ ! -e $SUBSETS ];then
    echo SUBSETS
    $MAIN --coverage --alignment-file $ALN > $SUBSETS || exit
fi

#rm -f $SUBTREES TRIPLETS $OUTDIR/parents $OUTDIR/build.tre $OUTDIR/strict.tre

echo SUBTREES
if [ ! -e $SUBTREES ];then
    $MAIN --display --parent-tree-file $TREE --subset-file $SUBSETS > $SUBTREES || exit
fi

echo TRIPLETS
if [ ! -e $TRIPLETS ];then
    $MAIN -t --subtree-file $SUBTREES > $TRIPLETS || exit
    #some giberish sometimes gets output to the stream, remove it
    perl -p -i -e 's/^.*1034h//g' $TRIPLETS
fi

echo PARENTS
if [ ! -e $OUTDIR/parentCount ];then
    $MAIN --count-parents --triplet-file $TRIPLETS > $OUTDIR/parentCount || exit
fi

echo GENERATE-PARENTS
PARENTS=$OUTDIR/allParents.tre
if [ ! -e $PARENTS ];then
    $MAIN --generate-parents --triplet-file $TRIPLETS > $PARENTS || exit
fi 


BUILD=$OUTDIR/build.tre 
if [ ! -e $BUILD ];then
    echo BUILD
    $MAIN --build --triplet-file $TRIPLETS > $BUILD || exit
fi


STRICT=$OUTDIR/strict.tre 
if [ ! -e $STRICT ];then
    echo STRICT
    $MAIN --strict --triplet-file $TRIPLETS > $STRICT || exit
fi


echo TERRACES
../scripts/terraphy.main.py --list-terraces --subset-file $SUBSETS --treefiles-to-assign $INPUTDIR/100pars.tre > $OUTDIR/terraceList

