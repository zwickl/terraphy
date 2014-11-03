#!/bin/bash

#input files
INPUTDIR=input/Bouchenak/
TREE=$INPUTDIR/bestTree.phy
ALN=$INPUTDIR/aln.exported.nex 

PERLSCRIPTDIR=perlScripts/

OUTPUTDIR=output/
mkdir -p $OUTPUTDIR

if [ -d "$PERLSCRIPTDIR" ];then
    $PERLSCRIPTDIR/parseNexusCharSetToDAV.pl -f $ALN > $OUTPUTDIR/perlSubsets || exit
    $PERLSCRIPTDIR/displaysub.pl $TREE $OUTPUTDIR/perlSubsets > $OUTPUTDIR/perlSubtrees || exit
    $PERLSCRIPTDIR/maketriplets.pl $OUTPUTDIR/perlSubtrees > $OUTPUTDIR/perlTriplets || exit
    $PERLSCRIPTDIR/countParents.pl $OUTPUTDIR/perlTriplets > $OUTPUTDIR/perlParents || exit
    #$PERLSCRIPTDIR/sameTerrace.pl $INPUTDIR/100pars.phy $OUTPUTDIR/perlSubsets $OUTPUTDIR/perlTerraceList || exit
fi

TERRAPHYDIR=../scripts/
MAIN=$TERRAPHYDIR/terraphy.main.py

#making the subsets currently takes a while in terraphy, so can copy over 
#the perl version if desired
#cp $OUTPUTDIR/perlSubsets $OUTPUTDIR/subsets
if [ ! -e $OUTPUTDIR/subsets ];then
    echo SUBSETS
    $MAIN --coverage -a $ALN > $OUTPUTDIR/subsets || exit
fi

rm -f $OUTPUTDIR/subtrees $OUTPUTDIR/triplets $OUTPUTDIR/parents $OUTPUTDIR/build.tre

echo SUBTREES
$MAIN --display --tree-files $TREE --subset-file $OUTPUTDIR/subsets > $OUTPUTDIR/subtrees || exit

echo TRIPLETS
$TERRAPHYDIR/maketriplets.py $OUTPUTDIR/subtrees > $OUTPUTDIR/triplets || exit
#grep -v "[ ].*[ ].*[ ].*" temp >> $OUTPUTDIR/triplets || exit
perl -p -i -e 's/^.*1034h//g' $OUTPUTDIR/triplets

echo PARENTS
$MAIN --parents --triplet-file $OUTPUTDIR/triplets > $OUTPUTDIR/parents || exit

echo BUILD
$MAIN --build --triplet-file $OUTPUTDIR/triplets > $OUTPUTDIR/build.tre || exit

echo TERRACES
../scripts/terraphy.main.py --list-terraces --subset-file $OUTPUTDIR/subsets --tree-files $INPUTDIR/100pars.tre > $OUTPUTDIR/terraceList

