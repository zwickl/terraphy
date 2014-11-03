#!/usr/bin/perl

# Generate the subtrees displayed by the given tree for the given taxa sets.
# Input: two files supplied on command line: first contains one newick tree, second contains k rows, where each row is a tab delim list of
# taxon names sampled from the original newick tree (ie, a subsets file).
# Output: k newick trees for the displayed subtrees

# USAGE: displaysub.pl treefile subsetsfile

###################

open FH, "<$ARGV[0]";
$nwk = <FH>;
close FH;

open FH, "<$ARGV[1]";
while (<FH>)
	{
	push @taxaSets, {map { $_ => 1 } split}; # array of the taxa sets (as hashes).
	}
close FH;



$numTaxaSets = @taxaSets;
for $ts (0..$numTaxaSets-1)
	{
	my $Tree = treeInit($nwk);
	@taxa = keys %{$taxaSets[$ts]};
	displayedSubtree($Tree,$taxaSets[$ts]); 	# pass ref to tree data structure and a ref to the taxa array
	writeNewick($Tree->{ROOT});
	}


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



