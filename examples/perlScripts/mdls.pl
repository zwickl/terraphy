#!/usr/bin/perl

# Solve the pairwise Maximum Definining Label Set problem.

# CURRENT IMPLEMENTATION SEEMS TO HAVE PROBLEMS UNLESS THE TWO TREES ARE ROOTED WITH A COMMON TAXON 
# AS THE SISTER GROUP OF EVERYTHING ELSE. THAT'S FINE, BECAUSE WE ALWAYS CAN ROOT WITH ONE OF THE
# COMMON TAXA AND REROOT LATER.

# USAGE: mdls.pl treefile

# Input is a file with two rooted binary trees, with at least some taxon overlap (2 taxa).

#############################

open FH, "<$ARGV[0]";
while (<>)
	{
	chomp;
	$s[$tree++] =<FH>;
	}
close FH;


$T1 = treeInit($s[0]);
$T2 = treeInit($s[1]);

#$T_OV = treeInit($s[2]);

$nwk_OV = makeIntersectTree($T1,$T2);
$T_OV = treeInit($nwk_OV);

print swriteNewick($T1->{ROOT}),"\n";
print swriteNewick($T2->{ROOT}),"\n";
print swriteNewick($T_OV->{ROOT}),"\n";


my %allTaxaH = allTaxa($T1,$T2);
$numUnion = scalar keys %allTaxaH;
print "Input trees processed. Sizes: $T1->{NLEAVES}, $T2->{NLEAVES}, overlap=$T_OV->{NLEAVES}, union=$numUnion\n";

buildTree($T1,$T2,$T_OV);

$numDiscarded = scalar @discarded;
print "Solution set has $maxTaxa leaves ($numDiscarded taxa discarded)\n";
print "Counts 00=$count00 10=$count10 01=$count01 11=$count11 SumCollapse=$sumWuda\n";

#foreach (@taxa)
#	{print "$_\n"} 

#print "@taxa\n";

sub allTaxa
{
my %allTaxaH;
foreach (@_)
	{
	my %leafH = %{$_->{LEAFH}};
	foreach my $t (keys %leafH)
		{$allTaxaH{$t}=1}
	}
return %allTaxaH
}

####################################################################
sub makeIntersectTree
{
my ($T1,$T2)=@_;

my %labels1=%{$T1->{LEAFH}};
my %labels2=%{$T2->{LEAFH}};

for $label1 (keys %labels1)
	{
	if (exists $labels2{$label1})
		{  $intersectTaxa{$label1}=1  }
	}


my %treeH;
my $root  = copyTree($T1->{ROOT});
$treeH{ROOT}=$root;
displayedSubtree(\%treeH,\%intersectTaxa);

return swriteNewick($treeH{ROOT});

#initIDTree(0,$root);
#initNLEAVES($root);
#my @tt = descLeafNodes($root);
#my %leafH = makeLeafHash(@tt);
#$treeH{LEAVES}=[@tt];
#$treeH{LEAFH}={%leafH};
#$treeH{NLEAVES}= $root->{NLEAVES};
#return \%treeH;
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

initIDTree(0,$root);
initNLEAVES($root);
my @tt = descLeafNodes($root);
my %leafH = makeLeafHash(@tt);
#foreach (keys %leafH) {print "$_\n"}
$treeH{ROOT}=$root;
$treeH{LEAVES}=[@tt];
$treeH{LEAFH}={%leafH};
$treeH{NLEAVES}= $root->{NLEAVES};
die ("Fatal error: tree is NOT binary\n") if (!isBinaryTree($root));
return \%treeH;
}

sub initIDTree
{
my ($startIX,$root)=@_;
$gIX=$startIX;
recurseIndexTree($root);
}
sub recurseIndexTree
{
my ($nodeRef)=@_;
$nodeRef->{ID}=$gIX++;
for $child (@{$nodeRef->{DESC}})
	{ recurseIndexTree($child); }
return ;
}
sub initNLEAVES
{
my ($nodeRef)=@_;
my $sum;
if (isLeaf($nodeRef)) {$sum=1} else {$sum=0}; 
for $child (@{$nodeRef->{DESC}})
	{ $sum += initNLEAVES($child); }
$nodeRef->{NLEAVES}=$sum;
return $sum;
}

sub invLCA
# returns two leaf nodes of given internal node, such that the node is the lca of the two leaves
{
my ($nd)=@_;
if (isLeaf($nd)) {return (undef,undef)}
$n1 = ${$nd->{DESC}}[0];
$n2 = ${$nd->{DESC}}[1];

while (!isLeaf($n1))
	{ $n1 = ${$n1->{DESC}}[0]; }
while (!isLeaf($n2))
	{ $n2 = ${$n2->{DESC}}[0]; }
return ($n1,$n2);

}

# Requires that there be a unique ID on each node of the tree. Finds path from node A to root, stores node IDs in a hash
# and then traverses back from node B until a match is found in the hash.

sub lca
{
my ($nodeA,$nodeB)=@_;

die "Missing args in lca function\n" if (!$nodeA || !$nodeB);

my @nodes=taxa_node_to_node($nodeA,0); # call to null will go all the way to the root
#print "ARGS: $nodeA->{ID}, $nodeB->{ID}\n";
#foreach (@nodes)
#	{print "lca array: $_->{ID}\n";}

my %nA = ();
foreach (@nodes)
	{ $nA{$_->{ID}}=1}


my $node = $nodeB;
while (! (exists $nA{$node->{ID}}))
	{
	if (isRoot($node)) {return $root};
	$node = $node->{ANC};
#print "$node->{ID}\n";
#die if ($count++ > 10);
	}
return $node;


}


####################################################################
sub taxa_node_to_node
# return a list of node refs from a shallow node to a deeper node, inclusive; stops if it gets to the root node and returns
# everything up to there including the root
{
my ($nodeShallow, $nodeDeep)=@_;
my ($node,@nodes);
my $node = $nodeShallow;
while ($node != $nodeDeep)
	{
	push @nodes, $node;
	if (isRoot($node)) {return @nodes};
	$node = $node->{ANC};
	}
return @nodes;
}

####################################################################
sub makeLeafHash
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

if (isLeaf($nodeRef))
	 {push @taxaT, $nodeRef;}
for $child (@{$nodeRef->{DESC}})
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
return {ID=>-1,NAME=>$name,DESC=>[],ANC=>-1,NLEAVES=>-1};
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
for $child (@{$nodeRef->{DESC}})
	{ if (!isBinaryTree($child)) {return 0}; }
return 1 ;
}
# **********************************************************

sub recursePrint
{
my ($nodeRef)=@_;
print "$nodeRef->{ID}: num leaves:$nodeRef->{NLEAVES}:taxon name = $nodeRef->{NAME}\n";
if (!isLeaf($nodeRef) && exists $nodeRef->{LCA1})
	{
	print "\tTree 1 invLCA: $nodeRef->{LCA1}{ID}\n";
	print "\tTree 2 invLCA: $nodeRef->{LCA2}{ID}\n";
	}
for $child (@{$nodeRef->{DESC}})
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
######################

sub offshootCount
	{
	# the sum of all the clade sizes of clades attached along the node path from A back rootward to B, not including A or B.
	# Also returns the taxa in question
	my ($a, $b)= @_;
	if (isRoot($a) || $a==$b) {return 0}
	die "OFFSHOOT ($a->{ID},$b->{ID}) not ancestor problem" if (!isA_anc_B($b,$a));
	my $sum=0;
	my @leafNames = ();
	my $n=$a->{ANC};
	#print "$n->{ID} $b->{ID}\n";
	my $prev=$a;
	#while ($n != $b)  
#...the following trap for isRoot keeps the code from hanging, but it shoould never get called...
	while ($n != $b && !isRoot($n))
		{
	#print "$n->{ID}\n";
		foreach $desc (@{$n->{DESC}})
			{
			next if ($desc == $prev); 	# skip the lineage back to the original node A of course
			$sum += $desc->{NLEAVES};
			push @leafNames, descLeafNames($desc);
	#print "$sum..$desc->{ID}\n";
			}
		$prev=$n;
		$n=$n->{ANC};
		}
	return ($sum,@leafNames);
	}
##########################
sub offshootCopy
	{
	# Adds offshoot subtrees along the path from srcA back rootward back to srcB, not including srcA or srcB, TO the edge
	# subtending node solA on the solution tree, which will have ancestor solB (to be determined).
	my ($solA,$srcA,$srcB)=@_;
	my $solB = $solA->{ANC};
	my $prevSol = $solA;
	my $nSrc = $srcA->{ANC};
	my $prevSrc = $srcA;
	while ($nSrc != $srcB  && !isRoot($n))
		{
		my $newSol = nodeNew("");
		intercalateNode($newSol,$prevSol,$solB);
		foreach my $descSrc (@{$nSrc->{DESC}})
			{
			next if ($descSrc == $prevSrc); 
			my $cpSrc = copyTree($descSrc);
			addChild($newSol,$cpSrc);
			}
		$prevSol=$newSol;
		$prevSrc=$prevSrc->{ANC};
		$nSrc=$nSrc->{ANC};
		}
	return
	}
##########################
sub intercalateNode
{
my ($node, $A,$ancA)=@_;
addChild($node,$A);
$A->{ANC}=$node;
$node->{ANC}=$ancA;
replaceChild($ancA,$A,$node);
}

sub replaceChild
{
my ($node,$replaceThis,$withThis)=@_;
for my $i (0..$#{$node->{DESC}})
	{
	if ($node->{DESC}[$i] == $replaceThis) {$node->{DESC}[$i]=$withThis}
	}
}
##########################
sub copyTree
	{
	my ($node)=@_;
	my $newNode = copyNode($node);
	foreach my $desc (@{$node->{DESC}})
		{
		my $newChild = copyTree($desc);
		addChild($newNode,$newChild);
		}
	return $newNode;
	}

sub copyNode
	# make new node with info copied from old: copies name and nleaves, but not ID number,anc or descs; 
	{
	my ($old)=@_;
	my $new = nodeNew($old->{NAME});
	$new->{NLEAVES}=$old->{NLEAVES};
	return $new;
	}

sub insertNode
	{
	my ($ndesc,$nanc)=@_;
	my $ninsert = nodeNew("");
	addChild($ninsert,$ndesc);
	addChild($nanc,$ninsert);
	return $ninsert;
	}

sub buildTree
	{
	my ($tree1,$tree2,$treeOV)=@_;
	my $root1=$tree1->{ROOT};
	my $root2=$tree2->{ROOT};
	my $rootOV=$treeOV->{ROOT};

	initializeLCAmaps($rootOV,$tree1,$tree2);
#	print "LCA Maps initialized...\n";
	$rootCopyTree = copyTree($rootOV);
	buildTreeRecurse($rootOV,$tree1,$tree2,$rootCopyTree);
	initIDTree(0,$rootCopyTree);
	initNLEAVES($rootCopyTree);
	my @tt = descLeafNodes($rootCopyTree);
	my %leafH = makeLeafHash(@tt);
#	recursePrint($rootCopyTree);
	writeNewick($rootCopyTree);
	}

sub initializeLCAmaps

	# at each internal node of the overlap tree, determine the corresponding LCA node on each of the input trees
	# and store the pair of lca nodes at this node of the overlap tree. 

		{
		my ($node,$T1,$T2)=@_;

#print "$node->{ID}\n";
		if (isLeaf($node)){return};
		my ($leaf1,$leaf2)=invLCA($node);
		my ($tax1)=$leaf1->{NAME};
		my ($tax2)=$leaf2->{NAME}; # these are two leaf labels on the OV tree that have this node as their LCA
#print "$tax1, $tax2\n";

	my $lca1 = lca($T1->{LEAFH}{$tax1}, $T1->{LEAFH}{$tax2}  );  

	my $lca2 = lca($T2->{LEAFH}{$tax1}, $T2->{LEAFH}{$tax2}  );  


	$node->{LCA1}=$lca1;
	$node->{LCA2}=$lca2;
	foreach my $desc (@{$node->{DESC}})
		{
		initializeLCAmaps($desc,$T1,$T2);
		}
	}

sub buildTreeRecurse
	{
	# For every edge of the OV tree, and its vertices, a,b, find the corresponding vertices in trees 1 and 2 and use that
	# to call the path comparison routine

	my ($node,$T1,$T2,$nodeCopyTree)=@_;
	if (!isRoot($node)) # for each edge (root edge doesn't exist)...
		{
		my ($n1,$n2);
		my $anc = $node->{ANC};
		my $nanc1 = $anc->{LCA1};
		my $nanc2 = $anc->{LCA2};
		if (isLeaf($node)) # if node is a leaf, we have to fetch node on input trees from leaf hash list, not lca maps
			{
			$n1 = $T1->{LEAFH}{$node->{NAME}};
			$n2 = $T2->{LEAFH}{$node->{NAME}};
			push @taxa,$node->{NAME};
			++$maxTaxa;			# we include the leaves of the overlap tree in the solution of course
			}
		else
			{
			$n1 = $node->{LCA1};
			$n2 = $node->{LCA2};
			}
#print "$node->{ID} ($node->{NAME}) ::  $n1->{ID}($n1->{NAME}):$nanc1->{ID} :: $n2->{ID}($n2->{NAME}):$nanc2->{ID}\n"; 

	my ($sum1,@leaves1) = offshootCount($n1,$nanc1);
	my ($sum2,@leaves2) = offshootCount($n2,$nanc2);
	#print "\n\nEdge below node $node->{ID}:\n";
	#print "Tree 1 number of taxa = $sum1 \n";
	#print "Tree 2 number of taxa = $sum2 \n";
	#print "Tree 1 number of taxa = $sum1 (@leaves1)\n";
	#print "Tree 2 number of taxa = $sum2 (@leaves2)\n";
	++$count00 if ($sum1==0 && $sum2==0);
	++$count10 if ($sum1>0 && $sum2==0);
	++$count01 if ($sum1==0 && $sum2>0);
	if ($sum1>0 && $sum2>0)
		{
		++$count11;
		$sumWuda += $sum1+$sum2; # this many taxa would be involved in a collapsed new node in a supertree;
		#...its a measure of how much better off we are by selectively removing some taxa.
		}

	if ($sum1 > $sum2) 
		{
		$maxTaxa += $sum1;
		push @taxa, @leaves1;
		push @discarded,@leaves2;
		offshootCopy ($nodeCopyTree,$n1,$nanc1);
		if ($sum2>0) {print_accepted($n1,$n2,1,\@leaves1,\@leaves2) }
		}
	if ($sum2 > $sum1) 
		{
		$maxTaxa += $sum2;
		push @taxa, @leaves2;
		push @discarded,@leaves1;
		offshootCopy ($nodeCopyTree,$n2,$nanc2);
		if ($sum1>0) {print_accepted($n1,$n2,2,\@leaves2,\@leaves1) }
		}
	if ($sum1 == $sum2 && $sum1>0)  # routine for ties: currently just uses tree 1!
		{
		$maxTaxa += $sum1;
		push @taxa, @leaves1;
		push @discarded,@leaves2;
		offshootCopy ($nodeCopyTree,$n1,$nanc1);
		print_accepted($n1,$n2,1,\@leaves1,\@leaves2);	
		}
	}
for my $ix (0..$#{$node->{DESC}})
	{
	buildTreeRecurse($node->{DESC}[$ix],$T1,$T2,$nodeCopyTree->{DESC}[$ix]);
	}
}
#############################################
sub print_accepted
{
my ($n1,$n2,$keepWhichTree,$leavesKeepRef,$leavesDiscardRef)=@_;

print "\nEdge between nodes $n1->{ID} and $n2->{ID}\n";
my $s1 = scalar @{$leavesKeepRef};
my $s2 = scalar @{$leavesDiscardRef};
print "\tKeeping $s1 taxa (from tree $keepWhichTree)...\t"; print "@{$leavesKeepRef}\n";
print "\tDiscarding $s2 taxa (from other tree)...\t"; print "@{$leavesDiscardRef}\n\n";
}
#############################################
sub isA_anc_B
{
my ($A,$B)=@_;
if ($A==$B) {return 0};
#print "In isanc: $A->{ID},$B->{ID}\n";
my ($n)=$B;
#print "ID:$n->{ID}\n";
while (!isRoot($n))
	{
	$n = $n->{ANC};
#print "ID:$n->{ID}\n";
	if ($n == $A) {return 1}
	}
return 0;
}
#############################################
sub writeNewick
{
my ($node)=@_;

if (isLeaf($node))
	{ print $node->{NAME} }
else
	{ 
	print "(";
	for my $ix (0..$#{$node->{DESC}})
		{
		writeNewick($node->{DESC}[$ix]);
		if ($ix < $#{$node->{DESC}}) { print "," };
		}

	print ")";
	}
if (isRoot($node)) {print ";\n"}
return;
}



###############################

sub displayedSubtree
{
my ($tRef,$taxaRef)=@_;

recursePrune($tRef->{ROOT},$taxaRef);

# SHOULD RECALCULATE THE NUMBER OF DESCENDANTS OF ALL NODES HERE...

return ;
}

sub recursePrune
{

# return 1 if the parent should KEEP this node as its child; 0 if the parent should delete this node as its child

my ($node,$taxaRef)=@_;


if (isLeaf($node))
	{
	if (exists $taxaRef->{$node->{NAME}}) 
		{return 1}
	else  
		{return 0}
	}
my @keptChildren=();
my $nDesc=scalar @{$node->{DESC}};
my $nKept=0;
foreach my $desc (@{$node->{DESC}})
	{
	if (recursePrune($desc,$taxaRef)==1)
		{
		push @keptChildren,$desc;
		++$nKept;
		}
	}
if ($nKept == $nDesc) # keep all descendants, no changes, so just return to parent signalling all's well
	{
	return 1;
	}
else
	{
	if ($nKept >= 2)
		{
		$node->{DESC} = [@keptChildren];
		return 1;
		}
	elsif ($nKept==1)  # we'll just make this node have the same information as its one descendant--effectively deleting the deg 1 node
		{
		my $onlyChild=$keptChildren[0];
		$node->{NAME}=$onlyChild->{NAME};
		$node->{DESC}=[@{$onlyChild->{DESC}}];
		return 1;
		}
	else 		# nKept = 0
		{
		return 0; # don't have to worry about the {DESC} data structure, we'll get rid of the node when we return anyway...
		}
	}
}

sub copyTree
        {
        my ($node)=@_;
        my $newNode = copyNode($node);
        foreach my $desc (@{$node->{DESC}})
                {
                my $newChild = copyTree($desc);
                addChild($newNode,$newChild);
                }
        return $newNode;
        }
sub copyNode
        # make new node with info copied from old: copies name,id, and nleaves, but not ID number,anc or descs;
        {
        my ($old)=@_;
        my $new = nodeNew($old->{NAME});
        $new->{NLEAVES}=$old->{NLEAVES};
        $new->{ID}=$old->{ID};
        return $new;
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

