package treesequal;

# Takes one file as a command line argument, which contains two newick strings
# ASSUMES TREES ARE ALL BINARY! Dies otherwise

# Uses the data structures of Day 1985 for O(n) strict consensus.

# The two tree structures, t1 and t2, are package globals and are accessed by various functions.

#use strict;

my ($t1,$t2);

sub treesequal
{

# reads two newick strings passed as args

$t1 = treeInit($_[0]);
min_maxLeafID($t1->{ROOT});
setupX($t1->{ROOT});

$t2 = treeInit($_[1]);
eval {
	min_maxLeafID_OTHER($t2->{ROOT}); # this will throw an exception if there is a label mismatch between trees! so that's good!
	if ($t1->{NLEAVES} != $t2->{NLEAVES} ) {die "Trees have different numbers of taxa\n"};
	compareTrees($t2->{ROOT});
};
deleteTree($t1);
deleteTree($t2);
if ($@)
	{ return "Trees are different: $@\n"; }
else
	{  return "Trees are the same\n"  }


}
####################################################################

sub compareTrees
{
my ($n)=@_;
if (isLeaf($n)) { return };

my $L2 = $n->{minDescID};
my $R2 = $n->{maxDescID};
die "Mismatched cluster size\n" if ($n->{NLEAVES}  !=  $R2 - $L2 + 1) ; # Exception! Mismatch based on cluster size

# OK, might be a matching cluster. If so, this pair, (L2,R2) will be stored in either row
# L2 or row R2 of the X matrix. Check for this and throw exception if not. See Day, p. 17.

die "Mismatched clade composition\n" 
if 	( !   
   		(
			( $XL[$L2] == $L2  && $XR[$L2] == $R2 )
					||
			( $XL[$R2] == $L2  && $XR[$R2] == $R2 )
		)
	);

for my $child (@{$n->{DESC}})
	{ compareTrees($child); }
return;
}


####################################################################
sub setupX
{
my ($n) = @_;
if (isLeaf($n)) {return}; # leave undefined for leaves
my $g;
if (isRoot($n))  # g(j) from eq(1) of Day 1985.
	{
	$g = $t1->{NLEAVES};
	}
else
	{
	my $thisIX = $n->{ID};
	my $nextIX = $thisIX + 1;  # f(j) + 1 
	my $nextNode = $t1->{NODES}[$nextIX];
	if (isLeaf($nextNode))
		{ $g = $n->{maxDescID}  } # select R of eq 1
	else
		{ $g = $n->{minDescID}  } # select R of eq 1
	}
$XL[$g]=$n->{minDescID};
$XR[$g]=$n->{maxDescID};
for my $child (@{$n->{DESC}})
	{ setupX($child); }
}

####################################################################
sub treeInit
{
my ($newick)=@_;
my ($root,%treeH);
$name = '[\w\d\_\.]+|\'.*?\''; # notice the non-greedy match needed in the '...' format
@tokens = ($newick=~/($name|\,|\(|\)|\;)/g);
$tokens_ix=0;
parseCheck();
$tokens_ix=0;
#print "@tokens\n";

my $tok = next_tok();
if ($tok =~ /\(/)
	{
	$root = make_group();
	}
else
	{die "First token is incorrect: $tok\n";}

initIDTree(1,1,$root);
initNLEAVES($root);
#min_maxLeafID($root);
undef @taxaAll; #used in descLeafNodes to return an array of all nodes indexed by the postorder node index $n->{ID}
my @tt = descLeafNodes($root);  # an array of all the leaf nodes in no particular order
my %leafH = makeLeafHash(@tt);  # keys are leaf names and values are nodes of the tree
#foreach (keys %leafH) {print "$_\n"}
$treeH{ROOT}=$root;
$treeH{LEAVES}=[@tt];
$treeH{NODES}=[@taxaAll];  # has size n+m+1 because we index from 1..n+m
$treeH{LEAFH}={%leafH};
$treeH{NLEAVES}= $root->{NLEAVES};



die ("Fatal error: tree is NOT binary\n") if (!isBinaryTree($root));
return \%treeH;
}

sub initIDTree # sets up postorder ids for all nodes [0,..,M-1], and leaf IDs on [0..N-1]. M nodes, N leaves
{
my ($startNodeIX,$startLeafIX,$root)=@_;
$gNodeIX=$startNodeIX;
$gLeafIX=$startLeafIX;
recurseIndexTree($root);
}
sub recurseIndexTree
{
my ($nodeRef)=@_;
if (isLeaf($nodeRef)) {$nodeRef->{leafID}=$gLeafIX++}
for my $child (@{$nodeRef->{DESC}})
	{ recurseIndexTree($child); }
$nodeRef->{ID}=$gNodeIX++; # these are the post order ids
return ;
}
sub initNLEAVES
{
my ($nodeRef)=@_;
my $sum;
if (isLeaf($nodeRef)) {$sum=1} else {$sum=0}; 
for my $child (@{$nodeRef->{DESC}})
	{ $sum += initNLEAVES($child); }
$nodeRef->{NLEAVES}=$sum;
return $sum;
}
####################################################################

sub min_maxLeafID # sets the min and max leaf Id for all the leaves descended from this internal node
{
my ($n)=@_;
if (isLeaf($n)) 
	{
	my $id = $n->{leafID};
	$n->{minDescID}=$id;
	$n->{maxDescID}=$id;
	return ($id,$id);
	}
my $tMin = +10000000; my $tMax = -1; 
for my $child (@{$n->{DESC}})
	{ 
	my ($min,$max) = min_maxLeafID($child); 
	if ($min < $tMin) {$tMin = $min} 
	if ($max > $tMax) {$tMax = $max} 
	}
$n->{minDescID}=$tMin;
$n->{maxDescID}=$tMax;
return ($tMin,$tMax);
}
####################################################################

sub min_maxLeafID_OTHER # sets the min and max leaf Id for tree but using TREE1's leaf ID's 
{
my ($n)=@_;
if (isLeaf($n))
        {
	my $thisName = $n->{NAME};		# name of this node on tree 2
	my $t1Node = $t1->{LEAFH}{$thisName};	# the corresponding node on tree 1
	if (!defined $t1Node)
		{die "Name not found on tree1 (trees different because of different label sets)\n"}
        my $id = $t1Node->{leafID};		#...and its leaf ID
        $n->{minDescID}=$id;
        $n->{maxDescID}=$id;
        return ($id,$id);
        }
my $tMin = +10000000; my $tMax = -1;
for my $child (@{$n->{DESC}})
        {
        my ($min,$max) = min_maxLeafID_OTHER($child);
        if ($min < $tMin) {$tMin = $min}
        if ($max > $tMax) {$tMax = $max}
        }
$n->{minDescID}=$tMin;
$n->{maxDescID}=$tMax;
return ($tMin,$tMax);
}



####################################################################
sub makeLeafHash  # a hash where keys are leaf names and values are nodes
{
my (%leafH)=();
foreach (@_)
	{$leafH{$_->{NAME}}=$_}
return %leafH;
}
####################################################################
sub descLeafNames
{
return map {$_->{NAME}} descLeafNodes(@_[0]);
}


sub descLeafNodes
{
my ($root)=@_;
undef @taxaT;
recursePush($root);
return @taxaT;
}

sub recursePush
{
my ($nodeRef)=@_;

$taxaAll[$nodeRef->{ID}] = $nodeRef;
if (isLeaf($nodeRef))
	 {push @taxaT, $nodeRef;}
for my $child (@{$nodeRef->{DESC}})
	{ recursePush($child); }
return ;
}



####################################################################
sub make_group
{
my $rootRef=nodeNew("");
while (my $tok = next_tok())
	{
	if ($tok =~ /$name/)
		{
		my $nodeRef = nodeNew($tok);
		addChild($rootRef,$nodeRef);
		}
	elsif ($tok =~ /\(/)
		{
		$nodeRef = make_group();
		addChild($rootRef,$nodeRef);
		}
	elsif ($tok =~ /\)/)
		{
		return $rootRef;
		}
	elsif ($tok =~ /,/)
		{
		next;
		}
	}
}


# **********************************************************

sub next_tok
{
if ($tokens_ix >= $#tokens) {return 0}
return ($tokens[$tokens_ix++])
}

# **********************************************************

sub nodeNew
{
my ($name)=@_;
return {ID=>-1,leafID=>-1,minDescID=>-1,maxDescID=>-1,NAME=>$name,DESC=>[],ANC=>-1,NLEAVES=>-1};
}
# **********************************************************
sub addChild 
{
my ($nodeRef,$childRef)=@_;
$childRef->{ANC}=$nodeRef;
push @{ ${$nodeRef}{DESC} },$childRef;
}


# **********************************************************
sub isBinaryTree
{
my ($nodeRef)=@_;
if (isLeaf($nodeRef)) {return 1};
if (scalar @{$nodeRef->{DESC}} != 2) {return 0};
for my $child (@{$nodeRef->{DESC}})
	{ if (!isBinaryTree($child)) {return 0}; }
return 1 ;
}
# **********************************************************

sub recursePrint
{
my ($n)=@_;
print "$n->{ID}: num leaves:$n->{NLEAVES}:taxon name = $n->{NAME}\tLeaf ID = $n->{leafID} MinDescID = $n->{minDescID} MaxDescId = $n->{maxDescID}\n";
for my $child (@{$n->{DESC}})
	{ recursePrint($child); }
return ;
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

#**********************************************************
sub parseCheck
{
my ($left,$right,$commas,$nTax);
for $tok (@tokens)
	{
	if ($tok =~ /\(/) {$left++};
	if ($tok =~ /\)/) {$right++};
	if ($tok =~ /$name/) {$nTax++};
	if ($tok =~ /\,/) {$commas++};
	}
#print "ntaxa=$nTax,left=$left, right=$right,commas=$commas\n";
die "Unmatched parens in newick string\n" if ($left != $right);
}
# **********************************************************
sub deleteTree
{
my ($TreeRef)=@_;

deleteTreeStructure($TreeRef->{ROOT});
undef @{$TreeRef->{LEAVES}};
undef @{$TreeRef->{NODES}};
undef %{$TreeRef->{LEAFH}};
undef $TreeRef->{ROOT};
undef $TreeRef;
}

sub deleteTreeStructure
{
my ($nodeRef)=@_;
for my $child (@{$nodeRef->{DESC}})
	{ deleteTreeStructure($child); }
undef $nodeRef->{DESC};
undef $nodeRef->{ID};
undef $nodeRef->{leafID};
undef $nodeRef->{minDescID};
undef $nodeRef->{maxDescID};
undef $nodeRef->{ANC};
undef $nodeRef->{NAME};
undef $nodeRef->{NLEAVES};
undef $nodeRef;
return ;
}
# **********************************************************


1;
