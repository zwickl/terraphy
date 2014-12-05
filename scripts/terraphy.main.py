#!/usr/bin/env python
import sys
import re
import os
import subprocess
import shlex
from argparse import ArgumentParser
from itertools import izip
from collections import Iterable

import cProfile, pstats, StringIO

sys.path.append("../")
from terraphy.triplets import *
from terraphy.dendroutils import compat_get_taxon_set, compat_encode_bipartitions, dendropy_read_treefile

#DENDROPY PACKAGE
try:
    from dendropy import TreeList, Tree, Node, DataSet, TaxonSet
    from dendropy.treecalc import symmetric_difference
    
    #for dendropy 4 compatability
    try:
        from dendropy.error import DataError
    except:
        from dendropy.utility.error import DataError


except ImportError:
    sys.exit('problem importing dendropy package - it is required')


#PYGRAPH PACKAGE
try:
    from pygraph.classes.graph import graph
    from pygraph.algorithms.accessibility import connected_components
except ImportError:
    sys.exit('problem importing pygraph package - it is required')


class MultiWriter(object):
    def __init__(self, streams):
        self.streams = streams if isinstance(streams, list) else [streams]

    def add_streams(self, more_streams):
        if isinstance(more_streams, list):
            self.streams.extend(more_streams)
        else:
            self.streams.append(more_streams)

    def write(self, output, more_streams=None):
        to_write = self.streams
        if more_streams:
            to_write += more_streams if isinstance(more_streams, list) else [more_streams]
        for stream in to_write:
            stream(output)


def profile_wrapper(func, profiler, *args, **kwargs):
    '''Just a wrapper to enable and disable profiling around
    a function call
    '''
    if profiler:
        profiler.enable()
    result = func(*args, **kwargs)
    if profiler:
        profiler.disable()
    return result


def calculate_triplets(intrees):
    triplets = []
    for tnum, tree in enumerate(intrees):
        treestr = tree.as_string(schema='newick', suppress_edge_lengths=True, suppress_rooting=True, suppress_internal_node_labels=True)
        #eliminate any singletons (tips without sisters), remove trailing whitespace or newlines
        treestr = re.sub('([A-Za-z_.-]+)', '"\\1"', treestr).rstrip()
        treestr = treestr.rstrip(';')
        #evaluate newick string as set of nested tuples
        triplets.extend(find_triplets_defining_edges_descending_from_node(eval(treestr)))
        sys.stderr.write('after tree %d: %d triplets\n' % (tnum, len(triplets)))

    #collect all labels included in any triplet, which should be the actual number 
    #included in the input trees
    all_taxa = set()
    for trip in triplets:
        for tax in trip:
            all_taxa.add(tax)

    return (all_taxa, triplets)


class CoverageMatrix(object):
    def __init__(self, rows=None):
        self.taxa = set()
        self.rows = {}
        self.columns = []

    def fill_from_subsets(self, subsets):
        self.columns = [set(col) for col in subsets]
        for col in self.columns:
            self.taxa |= col
        self.fill_rows()

    def __getitem__(self, index):
        if isinstance(index, str):
            if not self.rows:
                self.fill_rows()
            if index in self.rows:
                return self.rows[index]
            else:
                raise KeyError('Taxon %s not found in coverage matrix\n')
        else:
            assert(isinstance(index, int))
            if not self.columns:
                self.fill_columns()
            if len(self.columns) < index:
                return self.columns[index]

    def fill_columns(self, recalculate=False):
        if self.columns and not recalculate:
            sys.exit('columns already filled')
        elif not self.rows:
            sys.stderr.write('WARNING: filling matrix columns from empty matrix rows')
        for sub_num in xrange(len(self.rows[0])):
            self.columns.append([tax[sub_num] for tax in self.rows])

    def fill_rows(self, recalculate=False):
        if self.rows and not recalculate:
            sys.exit('taxa already filled')
        elif not self.columns:
            sys.stderr.write('WARNING: filling taxa from empty columns')
        self.rows = { tax:[] for tax in self.taxa }
        for col_num, col in enumerate(self.columns):
            col_set = set(col)
            for tax in self.taxa:
                if tax in col_set:
                    self.rows[tax].append(1)
                else:
                    self.rows[tax].append(0)
 

class IncompatibleTripletException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


VERBOSE=False
def build_or_strict_consensus(label_set, triplets, node, build):
    '''This function constructs the BUILD tree for given triplets or the strict consensus.
    The algorithms are essentially the same, the strict consensus just does a lot of
    extra work to see if an edge occurs in all trees before adding it.
    '''
    if not triplets:
        #No triplets, so no internal branches within clade.  So, add polytomy of leaves
        for label in label_set:
            node.add_child(Node(label=label))
    else:
        #This returns one component for each of the clades descending from this node
        #The members of the component indicate the labels in each clade
        #if there is only one component, some triplets are incompatible
        components = compute(label_set, triplets)
        if VERBOSE:
            print 'TO COMPUTE - build_or_strict'
            print label_set
            print triplets
            print components
            print '----'
        
        if len(components) > 1:
            for comp in components:
                if len(comp) == 1:
                    #if only one label in component, add leaf
                    #can't index sets, so need to use pop
                    node.add_child(Node(label=comp.pop()))
                elif len(comp) == 2:
                    #if two labels, add a cherry if the branch is in all trees, 
                    #otherwise add the two leaves to the current node
                    if VERBOSE:
                        print '\tTO is edge in all trees - build_or_strict = 2'
                        print '\t', comp
                        print '\t', triplets
                        print '\t----'
                    if build or is_edge_in_all_trees(comp, label_set, triplets):
                        new_node = node.add_child(Node())
                    else:
                        new_node = node
                    for el in comp:
                        new_node.add_child(Node(label=el))
                else:
                    #if > 2 labels in component, filter triplets to only 
                    #include those in which both ingroups and outgroup 
                    #are in the label set of the component, i.e. are in 
                    #the clade of interest
                    #then add an internal node for that clade and recurse
                    new_trip = winnow_triplets(comp, triplets)
                    if VERBOSE:
                        print '\tTO is edge in all trees - build_or_strict > 2'
                        print '\t', comp
                        print '\t', triplets
                        print '\t----'
                    if build or is_edge_in_all_trees(comp, label_set, triplets):
                        new_node = node.add_child(Node())
                    else:
                        new_node = node
                    build_or_strict_consensus(comp, new_trip, new_node, build)
        else:
            raise IncompatibleTripletException('Input is incompatible!')


def is_edge_in_all_trees(in_components, label_set, triplets):
    out_components = set(label_set)
    out_components -= in_components
    in_components = set(in_components)
    x = in_components.pop()
    for in_comp in in_components:
        for out_comp in out_components:
            for conflict in [set([(x, out_comp, in_comp)]), set([(in_comp, out_comp, x)])]:
                new_trip = triplets | conflict
                if VERBOSE:
                    print '\t\tTO are_triplets_compatible - is_edge_in_all_trees'
                    print '\t\tadded ', conflict
                    print '\t\t', label_set
                    print '\t\t', new_trip
                    print '\t\t----'
                if are_triplets_compatible(label_set, new_trip):
                    return False
    return True


def are_triplets_compatible(label_set, triplets):
    try:
        test_triplet_compatibility(label_set, triplets)
    except IncompatibleTripletException:
        if VERBOSE:
            print '\t\t\tincompat!'
        return False

    if VERBOSE:
        print '\t\t\tCOMPAT'
    return True


def test_triplet_compatibility(label_set, triplets):
    '''This is essentially the build algorithm, it just doesn't construct a tree
    Incompatibilty is indicated by raising an IncompatibleTripletException
    '''
    #print '', len(label_set), len(triplets)
    if not triplets:
        #No triplets, so no internal branches within clade.
        return
    else:
        #This returns one component for each of the clades descending from this node
        #The members of the component indicate the labels in each clade
        #if there is only one component, some triplets are incompatible
        if VERBOSE:
            print '\t\t\tTO COMPUTE - test_triplet_compatibility'
            print '\t\t\t', label_set
            print '\t\t\t', triplets
            print '\t\t\t----'
        components = compute(label_set, triplets)
        
        if len(components) > 1:
            for comp in components:
                if len(comp) == 1:
                    pass
                elif len(comp) == 2:
                    pass
                else:
                    #if > 2 labels in component, filter triplets to only 
                    #include those in which both ingroups and outgroup 
                    #are in the label set of the component, i.e. are in 
                    #the clade of interest
                    new_trip = winnow_triplets(comp, triplets)
                    test_triplet_compatibility(comp, new_trip)
        else:
            raise IncompatibleTripletException('blah')


def superb(label_set, triplets, num_parents):
    num_parents = 0
    if not triplets:
        num_parents = num_trees(len(label_set))

    else:
        components = compute(label_set, triplets)
        num_components = len(components)
        
        if num_components > 1:
            num_biparts = 2 ** (num_components - 1) - 1
            for i in xrange(1, num_biparts + 1):
                subset1, subset2 = create_bipartition(components, i)
                if len(subset1) <= 2:
                    q = 1
                else:
                    new_triplets = winnow_triplets(subset1, triplets)
                    q = superb(subset1, new_triplets, num_parents)

                if len(subset2) <= 2:
                    v = 1
                else:
                    new_triplets = winnow_triplets(subset2, triplets)
                    v = superb(subset2, new_triplets, num_parents)

                num_parents += q * v

        else:
            num_parents = 1

    return num_parents


def create_bipartition(components, nth):
    '''Take the components in the existing partition and divide them into two groups, 
    creating a biparition.  nth is an index that arbitarily generates one of the possible 
    groupings of the partition components into a bipartition, specifying by converting the 
    integer to a bit string
    NOTE: argument components MUST be ordered for this to work, so a list, not a set
    The returned subsets is a pair of sets.
    '''
    num_components = len(components)
    #bin will convert the integer to a binary string, starting with "0b". Slice that off.
    pattern = bin(nth)[2:]
    #leading zeros are trimmed, so add them if necessary to make the string the right length
    if len(pattern) < num_components:
        pattern = '0' * (num_components - len(pattern)) + bin(nth)[2:]
    subsets = (set( e for e in flattened_array_generator([components[el] for el in xrange(num_components) if pattern[el] == "0"])),
            set(e for e in flattened_array_generator([components[el] for el in xrange(num_components) if pattern[el] == "1"])))
    return subsets


def winnow_triplets(label_set, triplets):
    '''Return only those triplets for which all taxa appear in the label_set.
    label_set can be either a list or set, but this will be much faster as a set'''
    return { trip for trip in triplets if trip[0] in label_set and trip[1] in label_set and trip[2] in label_set }


def compute(label_set, triplets):
    '''Toward making a graph, add labels (graph nodes) as keys of dictionary, with 
    values being a set of other labels to which they are "connected", meaning that
    they appear together in a triplet as the ingroup.  Note that which ingroup label
    is noted as "connected" to the other is arbitrary.
    The triplet outgroups are ignored.
    label_set argument can be a list or set
    '''
    #labels are indicated as connected to themselves, which will be handy in my_connected_components
    connections = {lab:set([lab]) for lab in label_set}
    try:
        for in1, in2, out in triplets:
            connections[in1].add(in2)
    except ValueError:
        sys.exit('Could not unpack')
    except KeyError:
        raise

    #return pygraph_connected_components(connections)
    return my_connected_components(connections)


def my_connected_components(connections):
    '''This is my custom connected componenets implementation, which is much faster than the
    pygraph one, which probably has overhead due to its generality.
    Input is a dictionary of sets, with keys being taxon labels (nodes) and values being other 
    taxon labels (nodes) to which they are directly connected (i.e., are connected to by an edge).
    Assumes that each node has already been added as a value to it's own key.
    Note that the values may not be _ALL_ of the nodes to which the keys are connected.
    Each of the connections values are all of or part of a connected component.
    Return is a list of sets of connected componenets.
    '''
    assigned = set()
    components = []

    for star in connections.itervalues():
        #If any member of this partial component has already been put into a component (possibly
        #overlapping with multiple components) we'll need to collect anything that overlaps and
        #then remove the parts that were joined.  I don't see a faster way of doing this, except
        #maybe using frozensets for each component so that components can more easily be removed
        if star & assigned:
            tojoin = [ comp for comp in components if star & comp ]
            joined = star.union(*tojoin)
            for rem in tojoin:
                components.remove(rem)
            components.append(joined)
        else:
            components.append(star)
        assigned |= star

    return components


def pygraph_connected_components(connections):
    '''Using pygraph, make a graph (see compute function). Determine connected
    components, and then do a bit of munging to get into correct format.
    return value is a list of lists, with each inner list being one connected component.
    argument connections is a dictionary
    '''
    mygraph = graph()
    mygraph.add_nodes(connections.keys())
    for node, cons in connections.iteritems():
        for con in cons:
            mygraph.add_edge((node, con))

    connect_dict = connected_components(mygraph)

    #the connection dictionary looks like this, i.e. the nodes and the 
    #component number that they connect to 
    #{'A': 1, 'C': 1, 'B': 1, 'E': 2, 'D': 1, 'G': 2, 'F': 2, 'H': 2}
    #need to convert that into actual node sets for each componenet
    #it really doesn't seem like this should be necessary
    #comp = {num:[] for num in xrange(1, max(connect_dict.values()) + 1)}
    #for node, compnum in connect_dict.iteritems():
    #    comp[compnum].append(node)
    comp = {num:set() for num in xrange(1, max(connect_dict.values()) + 1)}
    for node, compnum in connect_dict.iteritems():
        comp[compnum].add(node)

    return comp.values()


def print_subsets(infile):
    dat = DataSet.get_from_stream(open(infile), schema='nexus')

    if not dat.char_matrices:
        sys.exit('no character data found')
    if len(dat.char_matrices) > 1:
        sys.exit('> 1 character matrix found')

    mat = dat.char_matrices[0]

    if not mat.character_subsets:
        sys.exit('no charsets found')

    taxon_labels = mat.taxon_set.labels()
    num_loci = len(mat.character_subsets.values())
    sub_matrices = []
    for val in mat.character_subsets.values():
        sub_matrices.append(mat.export_character_subset(val))
        
    presence_absence = dict.fromkeys(taxon_labels)

    for taxon in taxon_labels:
        presence_absence[taxon] = []
        for sub_mat in sub_matrices:
            if re.search('[acgt]', str(sub_mat[taxon]).lower()):
                presence_absence[taxon].append('1')
            else:
                presence_absence[taxon].append('0')

    #if options.verbose:
    #    for taxon in taxon_labels:
    #        print taxon, ''.join(presence_absence[taxon]), len([pa for pa in presence_absence[taxon] if pa == '1'])

    for locus in range(num_loci):
        print '\t'.join([tax for tax in presence_absence.keys() if presence_absence[tax][locus] == '1'])


def displayed_subtree(tree, labels, use_retain=False):
    #this is annoying, but Dendropy can consider the labels not matching depending on 
    #underscore vs. space issues
    #taxa = TaxonSet(tree.taxon_set)
    #newtree = Tree(tree, taxon_set=taxa)
    #DP3 vs. DP4 
    if hasattr(tree, 'taxon_namespace'):
        newtree = Tree(tree, taxon_namespace=tree.taxon_namespace)
    else:
        newtree = Tree(tree, taxon_set=tree.taxon_set)

    if isinstance(labels[0], str):
        if use_retain:
            newtree.retain_taxa_with_labels(labels)
        else:
            newtree.prune_taxa_with_labels(labels)
    else:
        if use_retain:
            newtree.retain_taxa(labels)
        else:
            newtree.prune_taxa(labels)

    compat_encode_bipartitions(newtree)
    return newtree


def print_displayed_subtrees(trees, subsets):
    sub_sets = [ set(subs) for subs in subsets ]
    all_taxa = set()
    for subs in sub_sets:
        all_taxa |= subs
    to_prune = [ list(all_taxa - subs) for subs in sub_sets ]

    for tnum, tree in enumerate(trees):
        tree.is_label_lookup_case_sensitive = True
        taxon_label_map = { taxon.label:taxon for taxon in compat_get_taxon_set(tree) }

        for setnum, (prune, retain) in enumerate(zip(to_prune, subsets)):
            if len(prune) < len(retain):
                sys.stderr.write('pruning tree %d to taxon set %d\n' % (tnum, setnum))
                newtree = displayed_subtree(tree, [ taxon_label_map[t] for t in prune ])
            else:
                sys.stderr.write('pruning tree %d to taxon set %d (using retain)\n' % (tnum, setnum))
                newtree = displayed_subtree(tree, [ taxon_label_map[t] for t in retain ], use_retain=True)
            if hasattr(newtree, 'as_newick_string'):
                newtreestr = newtree.as_newick_string() + ';\n'
            else:
                newtreestr = newtree.as_string(schema='newick', suppress_internal_node_labels=True, suppress_rooting=True)
            sys.stdout.write('%s' % newtreestr)


def same_tree(reference_tree, test_tree):
    '''This is adapted from false_positives_and_negatives() in dendropy treecalc module,
    and just bails when it finds the first different split.
    Should handle polytomies fine.
    '''

    ref_tset = compat_get_taxon_set(reference_tree)
    test_tset = compat_get_taxon_set(test_tree)
    if ref_tset is not test_tset:
        raise TypeError("Trees have different TaxonSet objects: %s vs. %s" \
                % (hex(id(ref_tset)), hex(id(test_tset))))

    compat_encode_bipartitions(reference_tree)
    compat_encode_bipartitions(test_tree)
    
    #seems like this set comparison should be faster, but not really
    if isinstance(reference_tree.split_edges, dict):
        if set(reference_tree.split_edges.keys()) != set(test_tree.split_edges.keys()):
            return False
    elif set(reference_tree.split_edges) != set(test_tree.split_edges):
        return False
    '''
    for split in reference_tree.split_edges:
        if split not in test_tree.split_edges:
            return False

    for split in test_tree.split_edges:
        if split not in reference_tree.split_edges:
            return False
    '''

    return True


def assign_to_terraces(trees, subsets):
    '''This is deprecated in favor of assign_to_terraces_using_hashes.'''
    #the first tree has to be its own terrace
    this_tree_subtrees = TreeList( [displayed_subtree(trees[0], subset) for subset in subsets] )
    terrace_subtree_list = [this_tree_subtrees]
    terrace_size = {0:1}
    
    #assign this same TaxonSet to all of the TreeLists created, otherwise bad things happen when trying to compare taxa and splits
    global_taxon_set = compat_get_taxon_set(this_tree_subtrees)
    
    #now start with the second tree
    for tree_num, tree in enumerate(trees[1:], 1):
        this_tree_subtrees = TreeList( [displayed_subtree(tree, subset) for subset in subsets], taxon_set=global_taxon_set )
        for terrace_num, terrace_subtrees in enumerate(terrace_subtree_list):
            terrace_match = True
            for sub1, sub2 in izip(terrace_subtrees, this_tree_subtrees):
                #if any subtree doesn't match, break to next terrace
                if not same_tree(sub1, sub2):
                    terrace_match = False
                    break
            #if all subtrees matched, terrace is same
            if terrace_match:
                print 'tree %d is on terrace %d' % (tree_num, terrace_num)
                terrace_size[terrace_num] += 1
                break
        else:
            #get here if went through all terraces with no match
            print 'tree %d is on a new terrace' % tree_num
            terrace_size[len(terrace_subtree_list)] = 1
            terrace_subtree_list.append(this_tree_subtrees)

    print terrace_size


def assign_to_terraces_using_hashes(trees, subsets):
    '''This is to be favored over assign_to_terraces.  The speed is about the same (somewhat surprisingly ...
    must be that the pruning function in displayed_subtree dominates the runtime), but since it 
    isn't actually storing any trees the memory usage doesn't grow as more terraces are found.
    '''
    #the first tree has to be its own terrace
    #this_tree_subtrees = TreeList( [displayed_subtree(trees[0], subset) for subset in subsets] )
    for tree in trees:
        compat_encode_bipartitions(tree)

    this_tree_subtrees = [ sum([ hash(n.edge.split_bitmask) for n in displayed_subtree(trees[0], subset).internal_nodes() ]) for subset in subsets]
    terrace_subtree_list = [this_tree_subtrees]
    terrace_size = {0:1}
    
    #now start with the second tree
    for tree_num, tree in enumerate(trees[1:], 1):
        this_tree_subtrees = [ sum([ hash(n.edge.split_bitmask) for n in displayed_subtree(tree, subset).internal_nodes() ]) for subset in subsets]
        for terrace_num, terrace_subtrees in enumerate(terrace_subtree_list):
            #if all subtrees matched, terrace is same
            if terrace_subtrees == this_tree_subtrees:
                print 'tree %d is on terrace %d' % (tree_num, terrace_num)
                terrace_size[terrace_num] += 1
                break
        else:
            #get here if went through all terraces with no match
            print 'tree %d is on a new terrace' % tree_num
            terrace_size[len(terrace_subtree_list)] = 1
            terrace_subtree_list.append(this_tree_subtrees)

    print terrace_size


def num_trees(taxa):
    '''The number of tree topologies for the given number of taxa'''
    trees = 0
    if taxa == 0:
        trees = 0
    elif taxa <= 2:
        trees = 1
    else:
        trees = 1
        for i in xrange(3, taxa + 1):
            trees *= ((2 * i) - 3)
    return trees


def draw_matrix_graphic(sorted_taxa, matrix, canvas):
    x_size = max(1, 1152 / len(matrix.rows))
    y_size = 50 / len(matrix.columns)
    x_loc, y_loc = 0, 0
    vertical = False
    if vertical:
        for tax, cov in matrix.rows.items():
            for cell in cov:
                if cell == 1:
                    canvas.create_rectangle(x_loc, y_loc, x_loc + x_size, y_loc + y_size, fill="blue", outline='blue')
                else:
                    canvas.create_rectangle(x_loc, y_loc, x_loc + x_size, y_loc + y_size, fill='white', outline='white')
                x_loc += x_size
            y_loc += y_size
            x_loc = 0
    else:
        for cov in matrix.columns:
            for tax in sorted_taxa:
                if tax in cov:
                    canvas.create_rectangle(x_loc, y_loc, x_loc + x_size, y_loc + y_size, fill="blue", outline='blue')
                else:
                    canvas.create_rectangle(x_loc, y_loc, x_loc + x_size, y_loc + y_size, fill='white', outline='white')
                x_loc += x_size
            y_loc += y_size
            x_loc = 0


########################################

parser = ArgumentParser(description='Perform various analyses related to phylogenetic terraces. Invoke script without arguments to start Tk GUI.')

in_group = parser.add_argument_group('Input Files')

in_group.add_argument('--alignment-file', default=None, help='nexus alignment including charsets to be used to determine character partition')

in_group.add_argument('--tree-files', nargs="*", default=None, help='nexus or newick tree file(s)')

in_group.add_argument('--subset-file', default=None, help='file with lines indicating sets of taxa represented in various partition subsets (created by --coverage preprocessing option)')

in_group.add_argument('--triplet-file', default=None, help='tab or space delimited triplet file, with ingroup ingroup outgroup (created by --triplets preprocessing option)')

preprocess = parser.add_argument_group('Preprosessing steps to perform on input files.  \nGeneral workflow would be --coverage, --display and --triplets, with each creating output consumed by following steps')

preprocess.add_argument('-c', '--coverage', action='store_true', default=False, help='compute the taxon coverage matrix (aka subsets file, requires --alignment-file)')

preprocess.add_argument('-d', '--display', action='store_true', default=False, help='print the subtrees displayed by the input tree with the input subsets (requires --subset-file and --tree-files)')

preprocess.add_argument('-t', '--triplets', action='store_true', default=False, help='Output arbitrary rooted taxon triples defining each edge in a set of treefiles (requires --tree-files')

analyses = parser.add_argument_group('Analyses to be performed on files created by preprocessing')

analyses.add_argument('-p', '--parents', action='store_true', default=False, help='compute the number of parent trees given a triplets file (requires --triplet-file')

analyses.add_argument('-b', '--build', action='store_true', default=False, help='compute the BUILD tree from a triplet file (requires --triplet-file)')

analyses.add_argument('-s', '--strict', action='store_true', default=False, help='compute a strict consensus tree from a triplet file (requires --triplet-file)')

analyses.add_argument('-l', '--list-terraces', action='store_true', default=False, help='take a set of trees and assign them to terraces (requires --subset-file and --tree-files)')

parser.add_argument('--profile', action='store_true', default=False, help='profile the given functionality')

writer = MultiWriter(sys.stdout.write)

#if no arguments are passed, try to start the tkinter gui
if len(sys.argv) == 1:
    try:
        from Tkinter import *
        from pygot.tkinterutils import *
        from ttk import *
    except ImportError:
        sys.stderr.write('%s\n' % parser.format_help())
        sys.stderr.write('\nUnable to import GUI componenets.  Use command line options.\n\n'.upper())
        sys.exit()

    tk_root = Tk()
    tk_gui = ArgparseGui(parser, tk_root, width=1280, height=720, output_frame=True, destroy_when_done=False, graphics_window=True)

    tk_root.mainloop()
    if tk_gui.cancelled:
        sys.exit('cancelled ...')
    options = parser.parse_args(tk_gui.make_commandline_list())
    output_result = tk_gui.output_result
    writer.add_streams(tk_gui.output_result)

else:
    tk_root = None
    options = parser.parse_args()
    output_result = sys.stdout.write

labels = []
triplets = set()
subsets = []
trees = None

prof = cProfile.Profile() if options.profile else None

if options.triplet_file:
    with open(options.triplet_file, 'rb') as intrips:
        #first line is a list of all labels, all following lines are triples of taxa
        for line in intrips:
            if not labels:
                labels = set(line.strip().split())
            else:
                trip = line.strip().split()
                #sort the ingroup taxa, to ensure no effectively identical triplets,
                #which would presumably come from different input trees
#                if trip[1] < trip[0]:
#                    trip[0], trip[1] = trip[1], trip[0]
                triplets.add(tuple(trip))

if options.subset_file:
    with open(options.subset_file, 'rb') as subs:
        for line in subs:
            subsets.append(line.strip().split())

    '''
    mat = CoverageMatrix()
    mat.fill_from_subsets(subsets)

    if tk_root:
        draw_matrix_graphic(mat.taxa, mat, tk_gui.graphics_canvas)
        
        def sort_and_redraw():
            taxa = list(mat.taxa)
            def cells(tax):
                return sum(mat.rows[tax])
            taxa.sort(key=cells, reverse=True)
            draw_matrix_graphic(taxa, mat, tk_gui.graphics_canvas)

        if hasattr(tk_gui, 'sort_button'):
            tk_gui.sort_button.config(command=sort_and_redraw)
        tk_root.mainloop()

    else:
        out_trans = ['-', 'X']
        for tax, cov in mat.rows.items():
            sys.stderr.write('%30s\t%s\n' % (tax, ''.join([out_trans[c] for c in cov])))
    '''


if options.tree_files:
    trees = dendropy_read_treefile(options.tree_files)

if options.triplets:
    all_taxa, triplets = profile_wrapper(calculate_triplets, prof, trees)
    sys.stdout.write('%s\n' % ' '.join([tax for tax in sorted(all_taxa)]))
    sys.stdout.write('%s\n' % '\n'.join([' '.join(trip) for trip in triplets]))

if options.build:
    if not options.triplet_file:
        sys.exit('triplet file (-t) must be supplied to make BUILD tree')
    
    build_tree = Tree()
    profile_wrapper(build_or_strict_consensus, prof, labels, triplets, build_tree.seed_node, build=True)
    writer.write('%s\n' % build_tree)
    
    #could automatically open in a viewer here
    #out_treefilename = 'build.tre'
    #with open(out_treefilename, 'w') as outtree:
    #    outtree.write('%s\n' % build_tree)
    #viewer_command =  shlex.split('java -jar "/Applications/FigTree1.3.1/FigTree v1.3.1.app/Contents/Resources/Java/figtree.jar"')
    #viewer_command.append(out_treefilename)
    #subprocess.call(viewer_command)
    
    if tk_root:
        tk_root.mainloop()

if options.strict:
    if not options.triplet_file:
        sys.exit('triplet file (-t) must be supplied to make strict consensus tree')
    
    strict_tree = Tree()
    profile_wrapper(build_or_strict_consensus, prof, labels, triplets, build_tree.seed_node, build=False)
    #strict(labels, triplets, strict_tree.seed_node)
    writer.write('%s\n' % strict_tree)
    
    #could automatically open in a viewer here
    #out_treefilename = 'strict.tre'
    #with open(out_treefilename, 'w') as outtree:
    #    outtree.write('%s\n' % strict_tree)
    #viewer_command =  shlex.split('java -jar "/Applications/FigTree1.3.1/FigTree v1.3.1.app/Contents/Resources/Java/figtree.jar"')
    #viewer_command.append(out_treefilename)
    #subprocess.call(viewer_command)
    
    if tk_root:
        tk_root.mainloop()

if options.parents:
    if not options.triplet_file:
        sys.exit('triplet file (-t) must be supplied to count parent tree')
    parents = profile_wrapper(superb, prof, labels, triplets, 0)
    output_result('%g\n' % parents)
    if tk_root:
        tk_root.mainloop()

if options.coverage:
    if not options.alignment_file:
        sys.exit('alignment file (-a) must be supplied to determine taxon coverage')
    profile_wrapper(print_subsets, prof, options.alignment_file)

if options.display:
    if not options.subset_file and options.tree_files:
        sys.exit('must specify both subset file (-s) and tree file (--tree-files) to print displayed subtrees')
    profile_wrapper(print_displayed_subtrees, prof, trees, subsets)

if options.list_terraces:
    if not options.subset_file and options.tree_files:
        sys.exit('must specify both subset file (-s) and tree file (--tree-files) to assign trees to terraces')
    profile_wrapper(assign_to_terraces, prof, trees, subsets)

if prof:
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(prof, stream=s).sort_stats(sortby)
    ps.print_stats()
    ps.print_callers()
    sys.stderr.write('%s\n' %s.getvalue())


