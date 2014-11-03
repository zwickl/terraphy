#!/usr/bin/perl

# For file of supplied trees, determine which trees are on distinct terraces
# Input: two files supplied on command line: first contains n>=2 newick trees, second contains k rows, 
# where each row is a tab delim list of taxon names sampled from the original newick tree.
# In other words, second file is a subsets file.

# YIKES. My first experience with PERL garbage collection, or should I say garbage strike.
# Had to write a formal destructor to get rid of trees and pruned parts of trees. Apparently,
# PERL's GC system has trouble with linked lists...

use lib '/usr/local/ActivePerl-5.14/site/lib/';
use lib 'perlScripts/';
die "USAGE: sameTerrace treefile subsetsfile outfile\n" if (3 != @ARGV);

use treesequal;

$tempfile = "sameisland.temp.$$"; # PERL process ID


open FH, "<$ARGV[1]";
while (<FH>)
	{
	push @taxaSets, {map { $_ => 1 } split}; # array of the taxa sets (as hashes).
	}
close FH;
$numTaxaSets = @taxaSets;

open FH1, ">$ARGV[2].$$";

open FH, "<$ARGV[0]";

$nwk1 = <FH>;
$ix1=0;
chomp $nwk1;

push @islands,$nwk1; # initialize the list of islands with the first tree

#print FH1 "Tree $ix1 is an island: adding to island list\n";
print "Tree $ix1 is an island: adding to island list\n";

while (<FH>)
	{
	++$ix1;
    #DJZ next tree string from file
	$nwk1=$_;
	chomp $nwk1;
	$numIslands = @islands;
	$exists_island=0;
    #DJZ iterate through the islands that have already been found
	for $ix2 (0..$#islands)
		{
        #DJZ tree string of one tree on the island?
		$nwk2=$islands[$ix2];
		$edge=1;
        #DJZ loop over taxon subsets
		for $ts (0..$numTaxaSets-1)
			{
            #DJZ take the new tree and one from the island being tested,
            #prune to this taxon subset and test for equality
			$Tree1 = treeInit($nwk1);
			$Tree2 = treeInit($nwk2);
			displayedSubtree($Tree1->{ROOT},$taxaSets[$ts]); 	
			displayedSubtree($Tree2->{ROOT},$taxaSets[$ts]); 
			$subnwk1[$ts]=swriteNewick($Tree1->{ROOT});
			$subnwk2[$ts]=swriteNewick($Tree2->{ROOT});
			#open FHO, ">$tempfile";
			#print FHO "$subnwk1[$ts]\n";
			#print FHO "$subnwk2[$ts]\n";
			#close FHO;
			#$s = "/home/sanderm/TREE_DEFINE/SAMEISLAND/treesequal.pl $tempfile";
			#$result = `$s` ;
			$result = treesequal::treesequal($subnwk1[$ts],$subnwk2[$ts]);
			deleteTree($Tree1); # important to do this before 'last' below
			deleteTree($Tree2);
            #DJZ result from treesequal comes back as a string?
			if ($result=~/different/)
				{
                #DJZ if any of the induced subtrees are different, trees aren't on same island
                #in that case, break and try the next island
				$edge = 0;
				last;
				}
			}
		if ($edge == 1) # this tree is found in an island already in the list
			{
            #DJZ for the island being tested all of the subtrees matched, so this tree must be on that island
            #keep track of that island number in saveIX2
			$exists_island=1; 
			$saveIX2=$ix2;
			last; 
			} # so we can bail on the rest of the island list
		}
	if ($exists_island == 0) # this is a new island,
		{ 
		push @islands, $nwk1; 
		#print FH1 "Tree $ix1 is a new island: adding to island list. Num islands =", $numIslands+1,"\n";
		print "Tree $ix1 is a new island: adding to island list. Num islands =", $numIslands+1,"\n";
		}
	else
		{
		#print FH1 "Tree $ix1 not a new island\n";
		print "Tree $ix1 not a new island: matches island $saveIX2 (0-off index)\n";
		}
	}

for $island (@islands)
	{
	print FH1 "$island\n";
	}

close FH1;


###############################

sub displayedSubtree
{
my ($root,$taxaRef)=@_;

recursePrune($root,$taxaRef);

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
	#***
	else
		{deleteTreeStructure($desc)}
	#***

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
sub deleteTree
{
my ($TreeRef)=@_;

deleteTreeStructure($TreeRef->{ROOT});
undef @{$TreeRef->{LEAVES}};
undef %{$TreeRef->{LEAFH}};
undef $TreeRef->{ROOT};
undef $TreeRef;
}

sub deleteTreeStructure
{
my ($nodeRef)=@_;
for $child (@{$nodeRef->{DESC}})
	{ deleteTreeStructure($child); }
undef $nodeRef->{DESC};
undef $nodeRef->{ID};
undef $nodeRef->{ANC};
undef $nodeRef->{NAME};
undef $nodeRef->{NLEAVES};
undef $nodeRef;
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
my ($left,$right,$nTax,$commas);
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
	# make new node with info copied from old: copies name,id, and nleaves, but not ID number,anc or descs; 
	{
	my ($old)=@_;
	my $new = nodeNew($old->{NAME});
	$new->{NLEAVES}=$old->{NLEAVES};
	$new->{ID}=$old->{ID};
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

