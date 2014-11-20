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
from terraphy.triplets import flattened_array_generator

#DENDROPY PACKAGE
try:
    from dendropy import TreeList, Tree, Node, DataSet, TaxonSet
    from dendropy.treecalc import symmetric_difference
    
    #for dendropy 4 compatability
    try:
        from dendropy.error import DataError
    except:
        from dendropy.utility.error import DataError

    #this deals with changes in DendroPy 4
    try:
        from dendropy.calculate import treesplit
    except ImportError:
        from dendropy import treesplit

except ImportError:
    sys.exit('problem importing dendropy package - it is required')


#PYGRAPH PACKAGE
try:
    from pygraph.classes.graph import graph
    from pygraph.algorithms.accessibility import connected_components
except ImportError:
    sys.exit('problem importing pygraph package - it is required')


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


def build(label_set, triplets, node):
    if not triplets:
        for label in label_set:
            node.add_child(Node(label=label))
    else:
        components = compute(label_set, triplets)
        
        if len(components) > 1:
            for comp in components:
                if len(comp) == 1:
                    #can't index sets, so need to use pop
                    node.add_child(Node(label=comp.pop()))
                elif len(comp) == 2:
                    new_node = node.add_child(Node())
                    for el in comp:
                        new_node.add_child(Node(label=el))
                else:
                    new_trip = winnow_triplets(comp, triplets)
                    new_node = node.add_child(Node())
                    build(comp, new_trip, new_node)
        else:
            sys.exit('triplets not compatible!')


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
    return [ trip for trip in triplets if trip[0] in label_set and trip[1] in label_set and trip[2] in label_set ]


def compute(label_set, triplets):
    '''Toward making a graph, add labels (graph nodes) as keys of dictionary, with 
    values being a set of other labels to which they are "connected", meaning that
    they appear together in a triplet as the ingroup.  Note that which ingroup label
    is noted as "connected" to the other is arbitrary.
    The triplet outgroups are ignored.
    label_set argument can be a list or set
    '''
    connections = {lab:set() for lab in label_set}
    try:
        for num, (in1, in2, out) in enumerate(triplets):
            connections[in1].add(in2)
    except ValueError:
        sys.exit('Could not unpack: %r' % triplets[0])

    #return pygraph_connected_components(connections)
    return my_connected_components(connections)


def my_connected_components(connections):
    assigned = set()
    components = []

    #connection_stars = [ set([node, *cons]) for node, cons in connections.iteritems() ]
    connection_stars = [ set([node] + list(cons)) for node, cons in connections.iteritems() ]

    #print 'stars', len(connection_stars)
    for star in connection_stars:
        if star & assigned:
            tojoin = []
            for comp in components:
                if star & comp:
                    tojoin.append(comp)
            #print 'comp, tojoin[0]', len(components), len(tojoin[0])
            joined = tojoin[0].union(star, *tojoin[1:])
            #print 'comp, tojoin[0]', len(components), len(tojoin[0])
            for rem in tojoin:
                components.remove(rem)
            #print 'comp', len(components)
            components.append(joined)
            #print 'comp', len(components)
            #components.append(tojoin[0])

        else:
            components.append(star)
        #print 'assigned', len(assigned)
        assigned |= star
        #print 'assigned', len(assigned)
        #print components
        #print 'assigned, all comp', len(assigned), len(components[0].union(*components[1:]))
        assert(assigned == components[0].union(*components[1:]))

    '''
    print 'len comp', len(components)
    for comp in components:
        for comp2 in components:
            if comp != comp2:
                #assert(not comp & comp2)
                if comp & comp2:
                    print comp & comp2
                    print sorted(list(comp))
                    print sorted(list(comp2))
    '''
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


def dendropy_read_treefile(filenames):
    intrees = TreeList()
    for tf in filenames:
        #try two input formats
        try:
            sys.stderr.write('Reading file %s in newick format ...\n' % tf)
            #intrees.extend(TreeList.get_from_path(tf, "newick"))
            intrees.extend(TreeList.get_from_path(tf, "newick", preserve_underscores=True))
        except DataError:
            sys.stderr.write('Reading file %s in nexus format ...\n' % tf)
            #intrees.extend(TreeList.get_from_path(tf, "nexus"))
            intrees.extend(TreeList.get_from_path(tf, "nexus", preserve_underscores=True))
        except ValueError:
            sys.exit('ValueError reading from file %s\n' % tf)
        except AttributeError:
            sys.exit('AttributeError reading from file %s\n' % tf)

    sys.stderr.write('read %d trees\n' % len(intrees))
    return intrees


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

    if hasattr(tree, 'encode_bipartitions'):
        tree.encode_bipartitions()
    else:
        treesplit.encode_splits(newtree)
    return newtree


def print_displayed_subtrees(trees, subsets):
    sub_sets = [ set(subs) for subs in subsets ]
    all_taxa = set()
    for subs in sub_sets:
        all_taxa |= subs
    to_prune = [ list(all_taxa - subs) for subs in sub_sets ]

    #labels = [ re.sub('_', ' ', label) for label in labels ]
    for tnum, tree in enumerate(trees):
        tree.is_label_lookup_case_sensitive = True
        #print type(tree.taxon_namespace[0].label)
        #taxon_label_map = { re.sub(' ', '_', taxon.label):taxon for taxon in tree.taxon_namespace }
        if hasattr(tree, 'taxon_namespace'):
            taxon_label_map = { taxon.label:taxon for taxon in tree.taxon_namespace }
        else:
            taxon_label_map = { taxon.label:taxon for taxon in tree.taxon_set }

        #print tree
        #print subsets
        #print to_prune
        for setnum, (prune, retain) in enumerate(zip(to_prune, subsets)):
            if len(prune) < len(retain):
                sys.stderr.write('pruning tree %d to taxon set %d\n' % (tnum, setnum))
                newtree = displayed_subtree(tree, [ taxon_label_map[t] for t in prune ])
                #newtree = displayed_subtree(tree, [ re.sub('_', ' ', taxon_label_map[t]) for t in prune ])
                #newtree = displayed_subtree(tree, prune)
            else:
                sys.stderr.write('pruning tree %d to taxon set %d (using retain)\n' % (tnum, setnum))
                newtree = displayed_subtree(tree, [ taxon_label_map[t] for t in retain ], use_retain=True)
                #newtree = displayed_subtree(tree, [ re.sub('_', ' ', taxon_label_map[t]) for t in retain ], use_retain=True)
                #newtree = displayed_subtree(tree, retain, use_retain=True)
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

    ref_tset = reference_tree.taxon_namespace if hasattr(reference_tree, 'taxon_namespace') else reference_tree.taxon_set
    test_tset = test_tree.taxon_namespace if hasattr(test_tree, 'taxon_namespace') else test_tree.taxon_set
    if ref_tset is not test_tset:
        raise TypeError("Trees have different TaxonSet objects: %s vs. %s" \
                % (hex(id(ref_tset)), hex(id(test_tset))))
    if hasattr(reference_tree, "encode_bipartitions"):
       reference_tree.encode_bipartitions()
    elif not hasattr(reference_tree, "split_edges") or not reference_tree.split_edges:
        treesplit.encode_splits(reference_tree)
    if hasattr(test_tree, "encode_bipartitions"):
       test_tree.encode_bipartitions()
    elif not hasattr(test_tree, "split_edges") or not test_tree.split_edges:
        treesplit.encode_splits(test_tree)
    
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
    global_taxon_set = this_tree_subtrees.taxon_namespace if hasattr(this_tree_subtrees, 'taxon_namespace') else this_tree_subtrees.taxon_set
    
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
        treesplit.encode_splits(tree)

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
    if (taxa == 0):
        trees = 0
    elif taxa <= 2:
        trees = 1
    else:
        trees = 1
        for i in xrange(3, taxa + 1):
            trees *= ((2 * i) -3)
    return trees


########################################

parser = ArgumentParser(description='Perform various analyses related to phylogenetic terraces. Invoke script without arguments to start Tk GUI.')

in_group = parser.add_argument_group('Input Files')

in_group.add_argument('-t', '--triplet-file', default=None, help='tab or space delimited triplet file, with ingroup ingroup outgroup')

in_group.add_argument('-s', '--subset-file', default=None, help='file with lines indicating sets of taxa represented in various partition subsets')

in_group.add_argument('-a', '--alignment-file', default=None, help='nexus alignment including charsets to be used to determine character partition')

in_group.add_argument('--tree-files', nargs="*", default=None, help='nexus or newick tree file(s)')

parser.add_argument('-b', '--build', action='store_true', default=False, help='compute the BUILD tree from a triplet file (requires --triplet-file)')

parser.add_argument('-p', '--parents', action='store_true', default=False, help='compute the number of parent trees given a triplets file (requires --triplet-file')

parser.add_argument('-c', '--coverage', action='store_true', default=False, help='compute the taxon coverage matrix (aka subsets file, requires --alignment-file)')

parser.add_argument('-d', '--display', action='store_true', default=False, help='print the subtrees displayed by the input tree with the input subsets (requires --subset-file and --tree-files)')

parser.add_argument('-l', '--list-terraces', action='store_true', default=False, help='take a set of trees and assign them to terraces (requires --subset-file and --tree-files)')

parser.add_argument('--profile', action='store_true', default=False, help='profile the given functionality')

#if no arguments are passed, try to start the tkinter gui
tk_root = None
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
    tk_gui = ArgparseGui(parser, tk_root, width=1152, height=720, output_frame=True, destroy_when_done=False)

    #Need to do this on OS X to bring window to front, otherwise root.lift() should work
    if 'darwin' in sys.platform.lower():
        try:
            #this can give odd non-critical error messages from the OS, so send stderr to devnull
            retcode = subprocess.call(shlex.split('''/usr/bin/osascript -e 'tell app "Finder" to set frontmost of process "Python" to true' '''), stderr=open(os.devnull, 'wb'))
        except:
            #didn't manage to get window to front, but don't worry about it
            pass
    else:
        tk_root.lift()
   
    tk_root.mainloop()
    if tk_gui.cancelled:
        sys.exit('cancelled ...')
    options = parser.parse_args(tk_gui.make_commandline_list())
    output_result = tk_gui.output_result
else:
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
                labels = line.strip().split()
            else:
                #triplets.append(line.strip().split())
                triplets.add(tuple(line.strip().split()))

if options.subset_file:
    with open(options.subset_file, 'rb') as subs:
        for line in subs:
            subsets.append(line.strip().split())

if options.tree_files:
    trees = dendropy_read_treefile(options.tree_files)

if options.build:
    if not options.triplet_file:
        sys.exit('triplet file (-t) must be supplied to make BUILD tree')
    
    build_tree = Tree()
    profile_wrapper(build, prof, labels, triplets, build_tree.seed_node)
    #build(labels, triplets, build_tree.seed_node)
    output_result('%s\n' % build_tree)
    
    #could automatically open in a viewer here
    #out_treefilename = 'build.tre'
    #with open(out_treefilename, 'w') as outtree:
    #    outtree.write('%s\n' % build_tree)
    #viewer_command =  shlex.split('java -jar "/Applications/FigTree1.3.1/FigTree v1.3.1.app/Contents/Resources/Java/figtree.jar"')
    #viewer_command.append(out_treefilename)
    #subprocess.call(viewer_command)
    
    if tk_root:
        tk_root.mainloop()

if options.parents:
    if not options.triplet_file:
        sys.exit('triplet file (-t) must be supplied to count parent tree')
    sys.stdout.write('%g\n' % profile_wrapper(superb, prof, labels, triplets, 0))

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


