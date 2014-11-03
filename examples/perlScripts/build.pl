#! /usr/bin/perl -w

# script to create the BUILD tree

# input: 
#	S: full label set
#	C: list of "constraints" (= rooted triplets)
# output:
#	nB: number of binary trees compatible with the constraints (binary parent trees)
#       treefile: BUILD tree, newick format
# v. 1.0 2011-03-12
#	Michelle McMahon, University of Arizona

use lib '/home/mcmahonm/bin/Graph-0.94/lib';
use lib '/usr/local/ActivePerl-5.14/site/lib/';
use strict; 
use warnings;
use Graph;
my $Usage = "

		./build.pl infile
		
		notes: 
		  1. infile must contain a list of taxa, space delimited, on one line
		     and a set of triples listed, one to a line, in the following format:
		         taxona taxonb taxonc
		     where taxonc is the sister to the other pair in the rooted triple
		  2. assumes (and doesn't check!) that input triplets are consistent
";

# ----- main -----
# -- input --
# get @S and @C from a file for now (can compute them later)
my $in = shift @ARGV;
die ("$Usage") if ((! $in) || (! -e $in));
my $h = 0;
my @S;
my @C = ();
open FHi, "<$in" or die "$!";
while (<FHi>)
  {
  next if (/^\s*$/);
  if (! $h) {@S = split; ++$h;}
  else 
    {
    my @temp = split; 
    die ("$!") if ((scalar @temp) != 3);
    push @C, [ @temp ];
    }
  }
close FHi;

# -- procedure --
#DJZ - the following two weren't used
#my $maxed = 0;
#my $totalTrees = 0;
my $Sref = \@S;
my $Cref = \@C;
my $gIX = 0;
my $rootRef = &nodeNew("root");
&BUILD($Sref, $Cref, $rootRef);
my $tree = &swriteNewick($rootRef);
print "$tree\n";

# ----- subroutine BUILD -----
sub BUILD
{
#DJZ These are the array of species names, the array of triples, and the root node
my ($Sref, $Cref, $Nref) = @_;
my @S = @{$Sref};
my @C = @{$Cref};
my $x = ${$Nref}{ID};

#DJZ if the triple array is empty, return star subtree
if (scalar @C == 0) 
  {
  #DJZ for each taxon name
  for my $S (@S) 
    {
    ++$gIX;
    my $newNodeRef = &nodeNew($S);
    &addChild($Nref, $newNodeRef);
    }
  } # end if C is empty
else
  {
  #DJZ get arrays of leaf labels that appear as ingroups in triplets together
  #i.e., given these triplets (with O=outgroup)
  #A B O
  #A C O
  #D E O
  #two sets in "partition" would be (A, B, C) and (D, E)
  my @pi_x = &COMPUTE($Sref, $Cref); # returns the partitions
  my $r = scalar @pi_x; # r = number of partitions
  
  #DJZ if there are at least 2 sets of triplets that don't intersect
  if ($r > 1)
    {
    for my $ii (0..($r-1)) # repeat across all partitions at this node
      {
      my $S1ref = $pi_x[$ii];
      my @S1 = @{$S1ref};
      #print "$ii @S1\n";
      #if the number of leaves in this set is 1
      if ((scalar @S1) == 1) 
        {
        #DJZ add the leaf
        ++$gIX;
        my $newNodeRef = &nodeNew($S1[0]);
        &addChild($Nref, $newNodeRef);
        }
      #DJZ if the number of leaves in this set is 2
      elsif ((scalar @S1) == 2) 
        {
        #DJZ add an internal node, followed by 2 desc
        ++$gIX;
        my $newNodeRef = &nodeNew("");
        &addChild($Nref, $newNodeRef);
        for my $S (@S1) 
          {
          ++$gIX;
          my $nextRef = &nodeNew("$S");
          &addChild($newNodeRef, $nextRef);
          }
        }
      #DJZ if there are more than 2 leaves in the set
      else
        {
        #DJZ get the set of triplets for which both ingroup leaves
        #appear in this label set
        my $newCRef = WINNOW_TRIPLETS($S1ref, $Cref);
        my @newC = @{$newCRef};
        ++$gIX;
        #DJZ add internal node
        my $newNodeRef = &nodeNew("");
        &addChild($Nref, $newNodeRef);
        #DJZ recursive call of this function using this set of labels as the 
        #new label list, the set of triplets that are fully included in that
        #label list as the triplet list, and the just added internal node
        #as the root node
        &BUILD($S1ref, $newCRef, $newNodeRef);
		}
      } # end for
    } # end if 
  else # r = 1
    {
    print "problem: triplets incompatible...\n";
    }
  } # end else (C not empty)

return;
} # end subroutine

# ----- subroutine WINNOW_TRIPLETS -----
sub WINNOW_TRIPLETS
{
my ($Xref, $Cref) = @_; # reference to label subset and full triple set
my %labH = ();
my $newCRef = [];
#DJZ dictionary of labels (used as a set, essentially)
for my $label (@{$Xref}) {$labH{$label} = 1;}
#DJZ loop over triplets
for my $trRef (@{$Cref})
  {
  my $fails = 0;
  #DJZ test whether all members of triplet appear in label subset
  for my $member (@{$trRef}) {if (! (exists $labH{$member})) {++$fails;} }
  next if ($fails); 
  push @{$newCRef}, $trRef;
  }
  #DJZ return only those triplets that are fully represented in the label set
return $newCRef;
}

# ----- subroutine COMPUTE -----
sub COMPUTE
{
my ($Sref, $Cref) = @_;
my @S = @{$Sref};
my @C = @{$Cref};
my $t; my $ii;
#DJZ dictionary of leaves that are connected to one another (=appear as
#the two ingroup taxa in some triplet) keys are leaves, values are other 
#leaves that they appear with in a triplet
my %connections;
for $t (@S) {$connections{$t}= [];}
for (@C) 
  {
    #DJZ for each ingroup pair in the triplets, use one as key and other
    #as value in dictionary, representing one edge. A single key may be
    #used multiple times, depending on the triples.  The outgroup isn't
    #used at all.
  my ($sis1, $sis2, $out) = @{$_};
  push @{$connections{$sis1}}, $sis2;
  }
my $href = \%connections;
my @subgraphs = &CONNECTED_COMPONENTS($href);
return (@subgraphs);
}

# ----- subroutine CONNECTED_COMPONENTS -----
sub CONNECTED_COMPONENTS
{
my ($href) = @_;
my %connections = %{$href};
my $g = Graph->new( undirected => 1 );
#DJZ key nodes are vertices of graph
$g->add_vertices(keys %connections);
for my $src ( keys %connections ) 
  {
#DJZ add edges to graph
  for my $tgt ( @{ $connections{$src} } ) {$g->add_edge($src, $tgt);}
  }
#return array of vertices that are connected (i.e., single linkage)
my @subgraphs = $g->connected_components;
#print "$subgraphs[0]\n";
return @subgraphs;
}

# **********************************************************
sub nodeNew
{
my ($name)=@_;
# ++$gIX;
return {ID=>$gIX,NAME=>$name,DESC=>[],ANC=>-1};
}
# **********************************************************
sub addChild 
{
my ($nodeRef,$childRef)=@_;
$childRef->{ANC}=$nodeRef;
push @{ ${$nodeRef}{DESC} },$childRef;
}
#############################################
sub swriteNewick
{
my ($node)=@_;

my $gs;

if (isLeaf($node))
         { $gs = $node->{NAME} }
else
         {
         $gs = "(";
         for my $ix (0..$#{$node->{DESC}})
                 {
                 $gs .= swriteNewick($node->{DESC}[$ix]);
                 if ($ix < $#{$node->{DESC}}) { $gs .= "," };
                 }

         $gs .= ")";
         }
if (isRoot($node)) {$gs .= ";"}
return $gs;
}

# **********************************************************

sub isLeaf
{
my ($nodeRef)=@_;
if (scalar @{$nodeRef->{DESC}} == 0)
	{return 1}
else
	{return 0}
}
sub isRoot
{
my ($nodeRef)=@_;
if ($nodeRef->{ANC} == -1)
	{return 1}
else
	{return 0}
}







