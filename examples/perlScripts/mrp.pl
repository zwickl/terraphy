#!/usr/bin/perl

# Makes an MRP matrix for a set of rooted trees

# USAGE: mrp.pl newickfile

# Takes one file as a command line argument, which contains newick strings
# Makes a ROOTED mrp matrix

##########################

open FH, "<$ARGV[0]";

while (<FH>)
	{
	chomp;
	$nwks[$ntree++] =$_;
	}
close FH;

$rootFlag=1;
foreach $nwk (@nwks)
	{
	$rootedTree = treeInit($nwk,$rootFlag);
	push @trees, $rootedTree;
	}
my %allTaxaH = allTaxa(@trees);
my @allTaxaA = sort keys %allTaxaH;
$numAllTaxa = scalar @allTaxaA;
my $ix=0;
foreach my $taxon (@allTaxaA)
	{
	$allTaxaIx{$taxon}=$ix++;
	}


$currentColumn=0;
foreach $tree (@trees)
	{
	makeMRP($tree);
	}
$ntax =$numAllTaxa;
$nchar=$currentColumn-1;
print "#nexus\nbegin data;\ndimensions ntax=$numAllTaxa nchar=$nchar;\nformat missing=? symbols=\"01\";\nmatrix\n";

for $ix (0..$ntax-1)
	{
	print "$allTaxaA[$ix]\t";
	for $iy (0..$nchar-1)
		{
		print "$mrp[$ix][$iy]";
		}
	print "\n";
	}
print ";\nend;\n";
#######################################
sub rootTree
{
my ($oldRoot)=@_;
my $og = nodeNew('__OG');
my $newRoot = nodeNew('');
addChild($newRoot,$og);
addChild($newRoot,$oldRoot);
return $newRoot;
}
####################################################################
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
sub makeMRP
{
my ($treeRef)=@_;
my $root = $treeRef->{ROOT};
my %leafNames = %{$treeRef->{LEAFH}};
@leafIndices= map { $allTaxaIx{$_}  } keys %leafNames; 
$numTaxaThisTree = scalar @leafIndices;
recurseMRP($root);
}


sub recurseMRP
	{
	my ($n) = @_;
	if (isLeaf($n)) {  return $n->{NAME}  };
	my @names = ();
	for my $child (@{$n->{DESC}})
		{ push @names, recurseMRP($child); }
	if (isRoot($n)) {return undef}
	my @cladeIndices = map { $allTaxaIx{$_}  } @names;	
	for my $ix (0..$numAllTaxa-1)
		{$mrp[$ix][$currentColumn]='?'}
	foreach my $ix (@leafIndices)
		{$mrp[$ix][$currentColumn]='0'}
	foreach my $ix (@cladeIndices)
		{$mrp[$ix][$currentColumn]='1'}
	++$currentColumn;
	return @names;
	}

####################################################################
sub treeInit
{
my ($newick,$rootFlag)=@_;
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
if ($rootFlag==1)
	{  $root = rootTree($root)  }
initIDTree(0,$root);
initNLEAVES($root);
my @tt = descLeafNodes($root);
my %leafH = makeLeafHash(@tt);
#foreach (keys %leafH) {print "$_\n"}
$treeH{ROOT}=$root;
$treeH{LEAVES}=[@tt];
$treeH{LEAFH}={%leafH};
$treeH{NLEAVES}= $root->{NLEAVES};
#die ("Fatal error: tree is NOT binary\n") if (!isBinaryTree($root));
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

sub oneDesc
# returns a single leaf descendant of a given internal node
{
my ($nd)=@_;
while (!isLeaf($nd))
        { $nd = ${$nd->{DESC}}[0]; }
return ($nd);

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
