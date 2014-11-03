#! /usr/bin/perl -w

# script to conduct the SUPERB algorithm of Constantinescu and Sankoff, 1995.
# lines with ## indicate text copied from appendix of paper. Contains a subroutine
# COMPUTE which is essentially the same as the BUILD algorithm of Aho et al..
# v 1.1 24 Feb 2011
# M. McMahon, Universty of Arizona

my $max = 1000000000;

# input: 
#	S: full label set
#	C: list of "constraints" (= rooted triplets)
# output:
#	nB: number of binary trees compatible with the constraints (binary parent trees)

use lib '/Users/sandermj/MyDOCUMENTS/PAPERS/TreeDefine.2010/CODE/lib/Graph-0.94/lib';
use lib '/usr/local/ActivePerl-5.14/site/lib/';
# use lib '/Users/mcmahonm/src/YAML-0.72/lib';
use strict; 
use warnings;
# use YAML;
use Graph;

my $Usage = "
		./countParents.pl infile
		
		notes: 
		  1. infile must contain a list of taxa, space delimited, on one line
		     and a set of triples listed, one to a line, in the following format:
		         taxona taxonb taxonc
		     where taxonc is the sister to the other pair in the rooted triple
		  2. assumes (and doesn't check!) that input triplets are consistent
		  3. max parents currently $max.
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
my $maxed = 0;
my $totalTrees = 0;
my $Sref = \@S;
my $Cref = \@C;
$totalTrees = &SUPERB($Sref, $Cref, $totalTrees);
if (! $maxed) {print "Total binary trees for input file $in: $totalTrees\n";}
else {print "Total binary trees for input file $in exceeds $max\n";}

# ----- subroutine SUPERB -----
sub SUPERB
{
my ($Sref, $Cref, $p) = @_;
my @S = @{$Sref};
my @C = @{$Cref};
my $s; my $q; my $v;

if (scalar @C == 0) 
  {
  #DJZ There are no triplets (constraints), so number of parents is total 
  #possible trees for that number of taxa
  my $ntaxa = scalar @S;
  $p = &NUMTREES($ntaxa);
  #print "no trips $p\n"
  } # end if C is empty
else
  {
      #print "normal\n";
  my @pi_x = &COMPUTE($Sref, $Cref); # returns the partitions
  my $r = scalar @pi_x; # r = number of partitions
  if ($r > 1)
    {
    $p = 0;
    $s = 2**($r-1)-1; # s = number of bipartitions
    #print "$r\t$s\n";

    for my $ii (1..$s) # repeat across all bipartitions at this node
      {
      my ($F1ref, $F2ref) = CREATE_BIPARTITION(\@pi_x, $ii);
      my @F1 = @{$F1ref}; 
      my @F2 = @{$F2ref};
      if ((scalar @F1) <= 2) {$q = 1;}
      else 
        {
        my $newCRef = WINNOW_TRIPLETS($F1ref, $Cref);
        my @newC = @{$newCRef};
        $q = SUPERB($F1ref, $newCRef, $p);
        }
      if ((scalar @F2) <= 2) {$v = 1;}
      else
        {
		my $newCRef = WINNOW_TRIPLETS($F2ref, $Cref);
        my @newC = @{$newCRef};
        $v = SUPERB($F2ref, $newCRef, $p);
        } # end else
      for my $t (1..$q)
        {
        for my $u (1..$v) {++$p; if ($p > $max) {$maxed=1; last;}} # ADDED MAX-TRAP HERE
        }
      } # end for
    } # end if 
  else # r = 1
    {
    $p = 1;
    }
  } # end else (C not empty)

return $p;
} # end subroutine

# ----- subroutine WINNOW_TRIPLETS -----
sub WINNOW_TRIPLETS
{
my ($Xref, $Cref) = @_; # reference to label subset and full triple set
my %labH = ();
my $newCRef = [];
for my $label (@{$Xref}) {$labH{$label} = 1;}
for my $trRef (@{$Cref})
  {
  my $fails = 0;
  for my $member (@{$trRef}) {if (! (exists $labH{$member})) {++$fails;} }
  next if ($fails); 
  push @{$newCRef}, $trRef;
  }
return $newCRef;
}

# ----- subroutine COMPUTE -----
sub COMPUTE
{
my ($Sref, $Cref) = @_;
my @S = @{$Sref};
my @C = @{$Cref};
my $t; my $ii;
my %connections;
for $t (@S) {$connections{$t}= [];}
for (@C) 
  {
  my ($sis1, $sis2, $out) = @{$_};
  push @{$connections{$sis1}}, $sis2;
  }
my $href = \%connections;
my @subgraphs = &CONNECTED_COMPONENTS($href);
return (@subgraphs);
}

# ----- subroutine NUMTREES -----
sub NUMTREES 
{
my ($n) = @_;
my $t = 0; my $ii;
die if ($n < 0);
if ($n == 0) {$t=0;}
elsif ($n <= 2) {$t=1;}
else 
  {
  $t=1;
  for $ii (3..$n) {$t = $t * (2*$ii -3);}
  }
return $t;
}

# ----- subroutine CONNECTED_COMPONENTS -----
sub CONNECTED_COMPONENTS
{
my ($href) = @_;
my %connections = %{$href};
my $g = Graph->new( undirected => 1 );
$g->add_vertices(keys %connections);
for my $src ( keys %connections ) 
  {
  for my $tgt ( @{ $connections{$src} } ) {$g->add_edge($src, $tgt);}
  }
my @subgraphs = $g->connected_components;
return @subgraphs;
}

# ----- subroutine CREATE_BIPARTITION -----
# DJZ - take the elements in the existing partition (objects) and
# divide them into two groups, creating a biparition.  Nth is an
# index that arbitarily generates one of the possible mappings of
# the partition into a bipartition, specifying by converting the 
# integer to a bit string
sub CREATE_BIPARTITION
{
my ($aref, $nth) = @_;
my @objects = @{$aref};
my $tot = scalar @objects;
my @partitions = ([ ], [ ]);
my @labels = ([ ], [ ]);
#DJZ This will be a bitstring representing the bipartition subset membership of each object
#in the initial partition
my $pattern = dec2bin($nth, $tot);
my @labels1 = (); my @labels2 = ();
#for each element of the original partition
for my $object (@objects)
  {
  #this gets the bipartition subset membership of each element of the initial partition 
  #from the next digit of the bitstring, and advances the string by replacing that digit
  #with ""
  my $part = substr($pattern,0,1,"");
  push @{$partitions[$part]}, @{$object};
  }
return (@partitions);
}

# ----- subroutine dec2bin -----
sub dec2bin 
{     
my ($dec, $tot) = @_;
my $str = unpack("B32", pack("N", $dec));
#DJZ grab only the last tot digits - everything
#before that will be zero
$str =~ s/^0+(\d{$tot})$/$1/;
return $str; 
}

