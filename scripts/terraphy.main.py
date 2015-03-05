#!/usr/bin/env python
import sys
import re
import os
import subprocess
import shlex
import threading
from argparse import ArgumentParser
from itertools import izip, combinations
from random import sample, random
from collections import Iterable
from copy import deepcopy

import cProfile, pstats, StringIO

sys.path.append("../")
from terraphy.triplets import *
from terraphy.coverage import *
from terraphy.dendroutils import compat_get_taxon_set, compat_encode_bipartitions, dendropy_read_treefile

#DENDROPY PACKAGE
try:
    from dendropy import TreeList, Tree, Node, DataSet, Taxon, TaxonSet
    #from dendropy.treecalc import symmetric_difference
    
    #for dendropy 4 compatability
    try:
        from dendropy.error import DataError
    except:
        from dendropy.utility.error import DataError

    try:
        from dendropy.treesim import uniform_pure_birth
    except:
        from dendropy.simulate.treesim import uniform_pure_birth_tree as uniform_pure_birth
    
    #this deals with changes in DendroPy 4
    try:
        from dendropy.calculate import treesplit
    except ImportError:
        from dendropy import treesplit

except ImportError as e:
    sys.exit('%s\nproblem importing dendropy package - it is required' % e)


#PYGRAPH PACKAGE
'''
#Implemented own connected components, so not currently using this
try:
    from pygraph.classes.graph import graph
    from pygraph.algorithms.accessibility import connected_components
except ImportError:
    sys.exit('problem importing pygraph package - it is required')
'''

#tree_viewer_command='java -jar "/Applications/FigTree1.3.1/FigTree v1.3.1.app/Contents/Resources/Java/figtree.jar"'
tree_viewer_command='java -jar "/Applications/FigTree1.3.1/lib/figtree.jar"'

class MultiWriter(object):
    '''A class that collects a number of output stream write methods and then calls all
    of them when desired.  Useful for example for sending output to Tkinter window
    and stdout without having to do cluter up code.
    Streams are actually the write methods of the given streams, i.e.
    stdout.write would be passed rather than stdout
    '''
    def __init__(self, streams):
        self.streams = streams if isinstance(streams, list) else [streams]

    def add_streams(self, more_streams):
        '''Append extra streams to always write to'''
        if isinstance(more_streams, list):
            self.streams.extend(more_streams)
        else:
            self.streams.append(more_streams)

    def write(self, output, more_streams=None):
        '''Call the collected write methods with the output string, plus any one time 
        streams passed in
        '''
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


def output_profile(prof):
    s = StringIO.StringIO()
    sortby = 'cumulative'
    ps = pstats.Stats(prof, stream=s).sort_stats(sortby)
    ps.print_stats()
    ps.print_callers()
    sys.stderr.write('%s\n' %s.getvalue())


def calculate_triplets(out, intrees, messages=sys.stderr):
    '''Choose a set of rooted triplets that define each edge of the input trees.  See the 
    terraphy.triplets find_triplets_defining_edges_descending_from_node function
    '''
    triplets = []
    all_taxa = set()
    if isinstance(intrees, str):
        intrees = dendropy_read_treefile([intrees], writer=messages)
    elif isinstance(intrees[0], str):
        intrees = dendropy_read_treefile(intrees, writer=messages)

    messages.write('generating triplets ..\n')
    for tnum, tree in enumerate(intrees):
        treestr = tree.as_string(schema='newick', suppress_edge_lengths=True, suppress_rooting=True, suppress_internal_node_labels=True)
        #eliminate any singletons (tips without sisters), remove trailing whitespace or newlines
        treestr = re.sub('([A-Za-z0-9_.-]+)', '"\\1"', treestr).rstrip()
        treestr = treestr.rstrip(';')
        #evaluate newick string as set of nested tuples
        tup_str = eval(treestr)
        all_taxa |= set( tax for tax in flattened_array_generator(tup_str, 9999999) )
        triplets.extend(find_triplets_defining_edges_descending_from_node(tup_str))
        messages.write('after tree %d: %d triplets\n' % (tnum, len(triplets)))
    
    messages.write('done.\n')
    
    if out:
        if isinstance(out, str):
            with open(out, 'w') as out_stream:
                out_stream.write('%s\n' % ' '.join([tax for tax in sorted(all_taxa)]))
                out_stream.write('%s\n' % '\n'.join([' '.join(trip) for trip in triplets]))
        else:
            out.write('%s\n' % ' '.join([tax for tax in sorted(all_taxa)]))
            out.write('%s\n' % '\n'.join([' '.join(trip) for trip in triplets]))
   
    #this is the legacy behavior, returning these instead of writing them in the func
    return (all_taxa, triplets)


def read_triplet_file(triplet_file, messages=sys.stderr):
    triplets = set()
    labels = []
    messages.write('Reading triplets from file... ')
    with open(triplet_file, 'rb') as intrips:
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
    messages.write('done.\n')
    return (labels, triplets)


class IncompatibleTripletException(Exception):
    '''An exception to allow escape from potentially deep recursion when looking for incompatible triplets'''
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


def debug_output(label_set, triplets, components, level=0):
    '''Specialized debug output for strict consensus computation'''
    indent = ''.join('\t' for l in xrange(level))
    if len(label_set) > 10:
        print indent, '%d labels, %d trips, %d comp' % (len(label_set or ''), len(triplets or ''), len(components or ''))
    else:
        if label_set:
            print indent, len(label_set), label_set
        if triplets:
            print indent, len(triplets), triplets
        if components:
            print indent, components
    print indent, '----'


def make_build_tree(out, triplet_file, messages=sys.stderr, verbose=False):
    tree = Tree()
    
    messages.write('Computing BUILD consensus tree...\n')
    label_set, triplets = read_triplet_file(triplet_file, messages=messages)
    
    build_or_strict_consensus(label_set, set(label_set), triplets, triplets, tree.seed_node, build=True, verbose=verbose)
   
    if out:
        if isinstance(out, str):
            with open(out, 'w') as out_stream:
                out_stream.write('%s;\n' % tree)
        else:
            out.write('%s;\n' % tree)
    
    messages.write('Done with BUILD consensus tree.\n')
    #this is the legacy behavior, returning instead of writing in the func
    return tree


def make_strict_tree(out, triplet_file, messages=sys.stderr, verbose=False):
    tree = Tree()
    
    messages.write('Computing strict consensus tree (this can take some time)...\n')
    label_set, triplets = read_triplet_file(triplet_file, messages=messages)
    
    build_or_strict_consensus(label_set, set(label_set), triplets, triplets, tree.seed_node, build=False, verbose=verbose)
   
    if out:
        if isinstance(out, str):
            with open(out, 'w') as out_stream:
                out_stream.write('%s;\n' % tree)
        else:
            out.write('%s;\n' % tree)
    
    messages.write('Done with strict consensus tree.\n')
    #this is the legacy behavior, returning instead of writing in the func
    return tree


def build_or_strict_consensus(label_set, full_label_set, triplets, all_triplets, node, build, precomp=None, verbose=False):
    '''This function constructs the BUILD tree for given triplets or the strict consensus.
    The algorithms are essentially the same, the strict consensus just does a lot of
    extra work to see if an edge occurs in all trees before adding it.
    '''

    if not build and not precomp:
        #if doing a strict consensus, make a dictionary of precomputed connected componenets, indexed by frozensets of 
        #label_sets. This will be used downstream in functions called by is_edge_in_all_trees
        precomp = {}
        compute_comp_dict(label_set, triplets, precomp)
    
    if not triplets:
        #No triplets, so no internal branches within clade.  So, add polytomy of leaves
        for label in label_set:
            node.add_child(Node(label=label))
    else:
        #This returns one component for each of the clades descending from this node
        #The members of the component indicate the labels in each clade
        #if there is only one component, some triplets are incompatible
        components = compute(label_set, triplets)
        if verbose:
            print 'AFTER COMPUTE - build_or_strict'
            debug_output(label_set, triplets, components)
        
        if len(components) > 1:
            for comp in components:
                if not comp:
                    print 'ZERO LENGTH COMPONENT'
                    sys.exit('ZERO LENGTH COMPONENT')
                if len(comp) == 1:
                    #if only one label in component, add leaf
                    #can't index sets, so need to use pop
                    node.add_child(Node(label=comp.pop()))
                    if verbose:
                        print '\tNEW SINGLETON'
                elif len(comp) == 2:
                    #if two labels, add a cherry if the branch is in all trees, 
                    #otherwise add the two leaves to the current node
                    if verbose:
                        print '\tTO is_edge_in_all_trees - build_or_strict = 2'
                        debug_output(comp, triplets, None, 1)
                    if build or is_edge_in_all_trees(comp, full_label_set, all_triplets, precomp=precomp, verbose=verbose):
                         if verbose:
                            print '\tNEW CHERRY'
                         new_node = node.add_child(Node())
                    else:
                        if verbose:
                            print '\tREJECTED', comp
                        new_node = node
                    for el in comp:
                        new_node.add_child(Node(label=el))
                else:
                    #if > 2 labels in component, filter triplets to only include those in which both ingroups and outgroup 
                    #are in the label set of the component, i.e. are in the clade of interest
                    #then if necessary add an internal node for that clade and recurse
                    new_trip = winnow_triplets(comp, triplets)
                    if verbose:
                        print '\tTO is_edge_in_all_trees - build_or_strict > 2'
                        debug_output(comp, triplets, None, 1)
                    if build or is_edge_in_all_trees(comp, full_label_set, all_triplets, precomp=precomp, verbose=verbose):
                         if verbose:
                            print '\tNEW CLADE'
                         new_node = node.add_child(Node())
                    else:
                        if verbose:
                            print '\tREJECTED', comp
                        new_node = node
                    build_or_strict_consensus(comp, full_label_set, new_trip, all_triplets, new_node, build, precomp=precomp, verbose=verbose)
        else:
            raise IncompatibleTripletException('Input is incompatible!')


def is_edge_in_all_trees(in_components, label_set, triplets, precomp=None, verbose=False):
    '''This is the heart of the Steel 1992 strict consensus algorithm.  To verify that an edge appears in all trees
    one needs to take the entire set of triplets, and iterate one by one over lots of other triplets that would conflict
    with the edge of interest.  The edge only appears in every tree if every one of those potentially conflicting triplets (PCTs)
    does conflict with the original set of triplets.  The procedure is:

    -Divide all labels into an ingroup set (those descending from the edge to be tested) and an outgroup set (all others)
    -Arbitrarily choose one of the ingroup labels to be a reference (x)
    -Iterate through all possible combinations of one of the remaining ingroup label (in_comp) and one outgroup label (out_comp)
    -For each combination of x, in_comp and out_comp, the real triplet resolution is obviously ((x, in_comp), out_comp), so the 
        PCTs are ((x, out_comp), in_comp) and ((in_comp, out_comp), x). Each is tested in turn by adding it to the set of real 
        triplets
    -To test for compatibility (which indicates that the edge of interest does NOT appear in all trees and therefore in the 
        strint consensus), a simplified version of the BUILD algorithm is used.  The set of triplets is recursively analyzed
        by calculating the connected components of a given set of triplets (clades descending from the current "node"), and
        recursing into each of those clades, after winnowing the triplets to only those contained in the clade.  The entire
        recursion starts with the entire set of triplets (at the root).
    -During the recursion, incompatibilty is indicated by only a single connected component of size > 2 existing at a given node

    Other notes:
    -If compatibilty is verified for any of the PCTs, the entire procedure is terminated, returning False from this function
    -During the testing of a particular PCT, once incompatibility is verfied at a given node the recursion is terminated
        for that PCT, moving on to the next one
    -This function only returns True if every single PCT is interated through and proven incompatible.

    There are a number of shortcuts that avoid full testing:
    -If a triplet exists on the same label set as a given PCT (i.e., ((x, in_comp), out_comp) is one of the real triplets)
        then recursion isn't necessary, incompatibilty is certain
    -During recursion, if the PCT is winnowed (i.e., not entirely contained within a clade), recursion into that clade
        isn't necessary (it can't conflict there)
    
    Optimizations:
    -In a naive implementation the most intensive computation comes from the computation of the connected components.
        However, they are determined directly from the set of triplets, and that set is exactly the same every time
        besides the addition of single PCTs.  So, it is possible to precompute the connected components induced by 
        a given label set (i.e., a given set of real triplets), and then just add the PCT to the computed components.
        That requires looking up the components induced by particular label sets, which is done with a dict indexed
        by frozensets of labels.
    '''

    out_components = set(label_set)
    out_components -= in_components
    in_components = set(in_components)
    x = in_components.pop()

    #this is the old method, precomputing only the root components
    #precomp = compute(label_set, triplets)
    
    #if one wasn't passed in, make a dictionary of precomputed connected componenets, indexed by frozensets of label_sets
    if not precomp:
        precomp = {}
        compute_comp_dict(label_set, triplets, precomp)

    for in_comp in in_components:
        for out_comp in out_components:
            #If a triplet for this exact taxon set exists, it will by definition conflict with the 
            #conflicting triplets added below.  So, don't bother doing the recursive checks.
            if set([(x, in_comp, out_comp), (in_comp, x, out_comp)]) & triplets:
                if verbose:
                    print '\t\tConflicting triplet already in triplet set, skipping recursion'
                next
            else:
                if verbose:
                    print '\t\tMust recurse to check for conflict'
                pass
            for conflict in [set([(x, out_comp, in_comp)]), set([(in_comp, out_comp, x)])]:
                new_trip = triplets | conflict
                if verbose:
                    print '\t\tTO are_triplets_compatible - is_edge_in_all_trees'
                    print '\t\tadded ', conflict
                    debug_output(label_set, new_trip, None, 2)
                
                #the new and old versions, with precomp of components at many nodes, or at only the root
                if isinstance(precomp, dict):
                    if are_triplets_compatible(label_set, new_trip, conflict, precomp=precomp, verbose=verbose):
                        return False
                else:
                    if are_triplets_compatible(label_set, new_trip, conflict, precomp=list(precomp), verbose=verbose):
                        return False
    return True


def are_triplets_compatible(label_set, triplets, conflict, precomp=None, verbose=False):
    '''Test compatibility of a set of triplets, starting a recursion into a function similar to BUILD.
    It is really more specialized than that, and assumes that a particular potentially conflicting triplet
    has already been added to the list of triplets, as well as passed as the conflict argument.  Used 
    in strict consensus computation. Incompatibility is indicated by catching an IncompatibleTripletException
    raised during recursion, compatibility by returning from recursion without an exception.
    May use precomputed connected components.  
    '''
    try:
        #call the root version, which may then recurse into normal test_triplet_compatibility, 
        #which doesn't use the precomputed components
        precomp_test_triplet_compatibility(label_set, triplets, conflict, precomp, verbose=verbose)
    except IncompatibleTripletException:
        if verbose:
            print 'INCOMPAT'
        return False

    if verbose:
        print '\t\t\tCOMPAT'
    return True


def test_triplet_compatibility(label_set, triplets, conflict, verbose=False, level=1):
    '''This is essentially the build algorithm, it just doesn't construct a tree
    Incompatibilty is indicated by raising an IncompatibleTripletException, which
    percolates upwards
    '''
    if not triplets:
        #No triplets, so no internal branches within clade.
        return
    else:
        #This returns one component for each of the clades descending from this node
        #The members of the component indicate the labels in each clade
        #if there is only one component, some triplets are incompatible
        components = compute(label_set, triplets)

        if verbose:
            indent = ''.join('\t' for l in xrange(level + 2))
            print indent, 'AFTER COMPUTE - test_triplet_compatibility'
            debug_output(label_set, triplets, components, level+2)

        if len(components) > 1:
            for comp in components:
                if not comp:
                    print 'ZERO LENGTH COMPONENT'
                    sys.exit('ZERO LENGTH COMPONENT')
                #these are the cases for a single tip or a cherry
                if len(comp) == 1:
                    if verbose:
                        print indent, 'COMP LEN 1 - PASS'
                    pass
                elif len(comp) == 2:
                    if verbose:
                        print indent, 'COMP LEN 2 - PASS'
                    pass
                else:
                    #if > 2 labels in component, filter triplets to only include those in which both ingroups and outgroup 
                    #are in the label set of the component, i.e. are in the clade of interest
                    
                    #doing the pre-winnowing test to see if the conflicting triplet will even be retained ends up being a bit
                    #faster in some cases here, although in cases where the total runtime is already short the overhead can make
                    #it slightly slower
                    #new_trip = winnow_triplets(comp, triplets)
                    #if conflict & new_trip:
                    if winnow_triplets(comp, conflict):
                        new_trip = winnow_triplets(comp, triplets)
                        if verbose:
                            print indent, 'COMP LEN %d - WINNOW AND RECURSE' % len(comp)
                        test_triplet_compatibility(comp, new_trip, conflict, verbose=verbose, level=level+1)
                    else:
                        if verbose:
                            print indent, 'COMP LEN %d - SKIP' % len(comp)
        else:
            if verbose:
                print ''.join('\t' for _ in  xrange(level+2)) ,
                #print 'INCOMPAT AT LEVEL %d' % level
                #print components
            raise IncompatibleTripletException('blah')


def precomp_test_triplet_compatibility(label_set, triplets, conflict, precomp=None, verbose=False, level=1):
    '''This is as test_triplet_compatibility, but allows for passing of precomputed components
    It is mainly a separate function to allow specific tracking during profiling
    '''
    if not triplets:
        #No triplets, so no internal branches within clade.
        return
    else:
        indent = ''.join('\t' for l in xrange(level + 2))
        #This returns one component for each of the clades descending from this node
        #The members of the component indicate the labels in each clade
        #if there is only one component, some triplets are incompatible
        if precomp:
            if isinstance(precomp, list):
                #This is the older way, where the only precomp was for the root
                if verbose:
                    print indent, 'Using list precomp'
                components = compute(set((list(conflict)[0][:2])), conflict, precomp)
            else:
                if frozenset(label_set) in precomp:
                    if verbose:
                        print indent, 'Using comp_dict'
                    components = compute(set((list(conflict)[0][:2])), conflict, precomp[frozenset(label_set)])
                else:
                    if verbose:
                        print indent, 'Not in comp_dict'
                    components = compute(label_set, triplets)
            
        else:
            components = compute(label_set, triplets)
        if verbose:
            print indent, 'AFTER COMPUTE - precomp_test_triplet_compatibility'
            if precomp:
                print indent, 'Used precomputed components'
            debug_output(label_set, triplets, components, level+2)

        if len(components) > 1:
            for comp in components:
                if not comp:
                    print 'ZERO LENGTH COMPONENT'
                    sys.exit('ZERO LENGTH COMPONENT')
                #these are the cases for a single tip or a cherry
                if len(comp) == 1:
                    if verbose:
                        print indent, 'COMP LEN 1 - PASS'
                    pass
                elif len(comp) == 2:
                    if verbose:
                        print indent, 'COMP LEN 2 - PASS'
                    pass
                else:
                    #if > 2 labels in component, filter triplets to only include those in which both ingroups and outgroup 
                    #are in the label set of the component, i.e. are in the clade of interest
                   
                    #doing the pre-winnowing test to see if the conflicting triplet will even be retained ends up being a bit
                    #faster in some cases here, although in cases where the total runtime is already short the overhead can make
                    #it slightly slower
                    #new_trip = winnow_triplets(comp, triplets)
                    #if conflict & new_trip:
                    if winnow_triplets(comp, conflict):
                        new_trip = winnow_triplets(comp, triplets)
                        if verbose:
                            print indent, 'COMP LEN %d - WINNOW AND RECURSE' % len(comp)
                        
                        if isinstance(precomp, dict) and frozenset(comp) in precomp:
                            precomp_test_triplet_compatibility(comp, new_trip, conflict, precomp=precomp, verbose=verbose, level=level+1)
                        else:
                            test_triplet_compatibility(comp, new_trip, conflict, verbose=verbose, level=level+1)
                    else:
                        if verbose:
                            print indent, 'COMP LEN %d - SKIP' % len(comp)
        else:
            if verbose:
                print ''.join('\t' for _ in  xrange(level+2)) ,
            raise IncompatibleTripletException('blah')


def compute_comp_dict(label_set, triplets, comp_dict, verbose=False, level=1):
    '''Precompute components for given label sets that can later be used as components to add conflicting triplets to during 
    strict consensus calculations'''
    
    if not triplets:
        return
    else:
        #This returns one component for each of the clades descending from this node
        #The members of the component indicate the labels in each clade
        #if there is only one component, some triplets are incompatible
        components = compute(label_set, triplets)
        comp_dict[frozenset(label_set)] = components
        if verbose:
            indent = ''.join('\t' for l in xrange(level))
            print indent, 'AFTER COMPUTE - compute_comp_dict'
            debug_output(label_set, triplets, components, level+2)

        if len(components) > 1:
            for comp in components:
                if not comp:
                    print 'ZERO LENGTH COMPONENT'
                    sys.exit('ZERO LENGTH COMPONENT')
                #these are the cases for a single tip or a cherry
                if len(comp) == 1:
                    pass
                elif len(comp) == 2:
                    pass
                else:
                    #if > 2 labels in component, filter triplets to only include those in which both ingroups 
                    #and outgroup are in the label set of the component, i.e. are in the clade of interest
                    new_trip = winnow_triplets(comp, triplets)
                    compute_comp_dict(comp, new_trip, comp_dict, verbose=verbose, level=level+1)
        else:
            if verbose:
                print '%d' % level ,
            raise IncompatibleTripletException('blah')


def count_trees_on_terrace(out, triplet_file, messages=sys.stderr):
    label_set, triplets = read_triplet_file(triplet_file, messages=messages)

    messages.write('Calculating parent trees on terrace...\n')
    num_par = superb(label_set, triplets)

    if out:
        if isinstance(out, str):
            with open(out, 'w') as out_stream:
                out_stream.write('%d trees on terrace\n' % num_par)
        else:
            out.write('%d trees on terrace\n' % num_par)
    
    return num_par


def superb(label_set, triplets):
    '''SUPERB algorithm of Constantinescu and Sankoff, 1995, to count number of parent trees
    compatible with given set of triplets'''
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
                    q = superb(subset1, new_triplets)

                if len(subset2) <= 2:
                    v = 1
                else:
                    new_triplets = winnow_triplets(subset2, triplets)
                    v = superb(subset2, new_triplets)

                num_parents += q * v

        else:
            num_parents = 1

    return num_parents


def create_bipartition(components, nth, as_list=False):
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
    if as_list:
        subsets = (list(subsets[0]), list(subsets[1]))
    return subsets


def winnow_triplets(label_set, triplets):
    '''Return a set containing only those triplets for which all taxa appear in the label_set.
    label_set can be either a list or set, but this will be much faster as a set'''
    
    #this is quite a bit faster than using "if label_set.issuperset(trip)"
    return { trip for trip in triplets if trip[0] in label_set and trip[1] in label_set and trip[2] in label_set }


def compute(label_set, triplets, precomp=None):
    '''Toward making a graph, add labels (graph nodes) as keys of dictionary, with 
    values being a set of other labels to which they are "connected", meaning that
    they appear together in a triplet as the ingroup.  Note that which ingroup label
    is noted as "connected" to the other is arbitrary.
    The triplet outgroups are ignored.
    label_set argument can be a list or set
    
    precomp is a list of (non-overlapping) sets of precomputed components, to which we want to add
    a potentially small number of new triplets. Label_set will only contain those labels appearing 
    as ingroups in the triplets passed in, not those in the precomputed components
    '''
    #labels are indicated as connected to themselves, which will be handy in my_connected_components
    connections = {lab:{lab} for lab in label_set}

    try:
        for in1, in2, out in triplets:
            connections[in1].add(in2)
    except ValueError:
        sys.exit('Could not unpack')
    except KeyError:
        raise

    #return pygraph_connected_components(connections)
    return my_connected_components(connections, precomp)


def my_connected_components(connections, precomp=None):
    '''This is my custom connected componenets implementation, which is much faster than the
    pygraph one, which probably has overhead due to its generality.
    Input is a dictionary of sets, with keys being taxon labels (nodes) and values being other 
    taxon labels (nodes) to which they are directly connected (i.e., are connected to by an edge).
    Assumes that each node has already been added as a value to it's own key.
    Note that the values may not be _ALL_ of the nodes to which the keys are connected.
    Each of the connections values are all of or part of a connected component.
    Return is a list of sets of connected componenets.

    A precomputed (partial) set of connected components can be passed in, and any additional
    connections are added to it.  A complication is that if the additional connections do 
    alter the precomp connected components (by joining some), then the actual precomp list
    would be changed.  So, in that case the precomp must be deepcopied and altered before
    being returned.
    '''

    if precomp:
        #components have been precomputed, but we want to add some connections
        components = precomp
        assigned = precomp[0].union(*precomp[1:])
        dcopied = False
    else:
        components = []
        assigned = set()
    
    if False:
        #this ends up being slower than the old version below, as do other tweaks I tried
        for star in connections.itervalues():
            if star & assigned:
                to_remove = []
                for comp in components:
                    if comp & star:
                        star |= comp
                        to_remove.append(comp)
                for comp in to_remove:
                    components.remove(comp)
                components.append(star)
            else:
                components.append(star)
            assigned |= star
    else:
        for star in connections.itervalues():
            #If any member of this partial component has already been put into a component (possibly
            #overlapping with multiple components) we'll need to collect anything that overlaps and
            #then remove the parts that were joined from the component list.  I don't see a faster 
            #way of doing this, except maybe using frozensets for each component so that components 
            #can more easily be removed
            if star & assigned:
                tojoin = [ comp for comp in components if star & comp ]
                joined = star.union(*tojoin)
                
                #if this isn't true, then the existing components won't be changed by including this star
                if len(tojoin) > 1 or joined != tojoin[0]:
                    if precomp and not dcopied:
                        #In this case the precomputed components (list of sets) will have to be 
                        #altered, which would change their values in the comp dictionary and screw 
                        #up their use in the future. So, copy to a new list that will be returned
                        #print "COPIED"
                        components = deepcopy(components)
                        dcopied = True
                    #seems like the listcomp might be faster here, but it isn't
                    #components = [ comp for comp in components if comp not in tojoin ]
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


def print_subsets(out, infile, messages=sys.stderr):
    messages.write('Reading alignment file ...\n')
    if isinstance(infile, list):
        infile = infile[0]
    dat = DataSet.get_from_stream(open(infile), schema='nexus')

    if not dat.char_matrices:
        sys.exit('no character data found')
    if len(dat.char_matrices) > 1:
        sys.exit('> 1 character matrix found')

    mat = dat.char_matrices[0]

    if not mat.character_subsets:
        sys.exit('no charsets found')

    messages.write('Calculating coverage matrix ... ')
    char_subsets = mat.character_subsets.values()
    
    presence_absence = {}

    if isinstance(out, str):
        out_stream = open(out, 'w')
    else:
        out_stream = out
    for taxon in mat.taxon_set:
        presence_absence[taxon] = []
        seq = mat.taxon_seq_map[taxon]

        for subset_indeces in char_subsets:
            for index in subset_indeces:
                if str(seq[index]) in 'acgtACGT':
                    presence_absence[taxon].append(True)
                    break
            else:
                presence_absence[taxon].append(False)

    for locus in range(len(char_subsets)):
        #print '\t'.join([re.sub(' ', '_', tax.label) for tax in presence_absence.iterkeys() if presence_absence[tax][locus]])
        out_stream.write('%s\n' % '\t'.join([re.sub(' ', '_', tax.label) for tax in presence_absence.iterkeys() if presence_absence[tax][locus]]))
    
    if isinstance(out, str):
        out_stream.close()
    
    messages.write('done.\n')


def read_subset_file(subset_file):
    subsets = []
    with open(subset_file, 'rb') as subs:
        for line in subs:
            sub = line.strip().split()
            #ignore blank lines
            if sub:
                subsets.append(sub)
    return subsets
 

def displayed_subtree(tree, labels, use_retain=False):
    #this is annoying, but Dendropy can consider the labels not matching depending on 
    #underscore vs. space issues
    #Realized that using retain_taxa vs prune_taxa doesn't make any difference - dendropy
    #used prune behind the scenes regardless
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

    #not sure why I had delete_outdegree_one here, but doesn't work with DP4
    #compat_encode_bipartitions(newtree, delete_outdegree_one=False)
    compat_encode_bipartitions(newtree)
    return newtree


def print_displayed_subtrees(out, trees, subsets, messages=sys.stderr):
    '''To be more flexible with callbacks, etc., allowing input and output to either be
    filenames (if strings) or objects.
    '''
    messages.write('Computing subtrees ...\n')
    if isinstance(out, str):
        out_stream = open(out, 'w')
    else:
        out_stream = out

    if isinstance(trees, str):
        trees = dendropy_read_treefile([trees], writer=messages)
    elif isinstance(trees[0], str):
        trees = dendropy_read_treefile(trees, writer=messages)

    if isinstance(subsets, str):
        messages.write('Reading subset file ... ')
        subsets = read_subset_file(subsets)
        messages.write('done.\n')
    elif isinstance(subsets[0], str):
        messages.write('Reading subset file ... ')
        subsets = read_subset_file(subsets[0])
        messages.write('done.\n')
    
    sub_sets = [ set(subs) for subs in subsets ]
    all_taxa = set()
    for subs in sub_sets:
        all_taxa |= subs
    
    for tnum, tree in enumerate(trees):
        tree.is_label_lookup_case_sensitive = True
        taxon_label_map = { taxon.label:taxon for taxon in compat_get_taxon_set(tree) }

        for setnum, retain in enumerate(subsets):
            messages.write('pruning tree %d to taxon set %d\n' % (tnum, setnum))
            newtree = displayed_subtree(tree, [ taxon_label_map[t] for t in retain ], use_retain=True)
            newtreestr = newtree.as_string(schema='newick', suppress_internal_node_labels=True, suppress_rooting=True)
            out_stream.write('%s' % newtreestr)
       
    messages.write('done.\n')


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


def assign_to_terraces(out, treefiles, subset_file, messages=sys.stderr):
    '''This is deprecated in favor of assign_to_terraces_using_hashes.'''
    
    messages.write('Assigning trees to terraces...\n')
 
    if isinstance(treefiles, str):
        trees = dendropy_read_treefile([treefiles], writer=messages)
    elif isinstance(treefiles[0], str):
        trees = dendropy_read_treefile(treefiles, writer=messages)
   
    subsets = read_subset_file(subset_file)
    
    #the first tree has to be its own terrace
    this_tree_subtrees = TreeList( [displayed_subtree(trees[0], subset) for subset in subsets] )
    terrace_subtree_list = [this_tree_subtrees]
    terrace_size = {0:1}
   
    if out:
        if isinstance(out, str):
            out_stream = open(out, 'w')
        else:
            out_stream = out

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
                out_stream.write('tree %d is on terrace %d\n' % (tree_num, terrace_num))
                terrace_size[terrace_num] += 1
                break
        else:
            #get here if went through all terraces with no match
            out_stream.write('tree %d is on a new terrace\n' % tree_num)
            terrace_size[len(terrace_subtree_list)] = 1
            terrace_subtree_list.append(this_tree_subtrees)

    messages.write(' done.')
    out_stream.write('%d terraces found' % len(terrace_size))


def assign_to_terraces_using_hashes(out, treefiles, subset_file, messages=sys.stderr):
    '''This is to be favored over assign_to_terraces.  The speed is about the same (somewhat surprisingly ...
    must be that the pruning function in displayed_subtree dominates the runtime), but since it 
    isn't actually storing any trees the memory usage doesn't grow as more terraces are found.
    '''
    messages.write('Assigning trees to terraces...\n')
    
    if isinstance(treefiles, str):
        trees = dendropy_read_treefile([treefiles], writer=messages)
    elif isinstance(trees[0], str):
        trees = dendropy_read_treefile(treefiles, writer=messages)
   
    for tree in trees:
        compat_encode_bipartitions(tree)

    subsets = read_subset_file(subset_file)
    
    if out:
        if isinstance(out, str):
            out_stream = open(out, 'w')
        else:
            out_stream = out
 
    #the first tree has to be its own terrace
    #this_tree_subtrees = TreeList( [displayed_subtree(trees[0], subset) for subset in subsets] )
    this_tree_subtrees = [ sum([ hash(n.edge.split_bitmask) for n in displayed_subtree(trees[0], subset).internal_nodes() ]) for subset in subsets]
    terrace_subtree_list = [this_tree_subtrees]
    terrace_size = {0:1}
    
    #now start with the second tree
    for tree_num, tree in enumerate(trees[1:], 1):
        this_tree_subtrees = [ sum([ hash(n.edge.split_bitmask) for n in displayed_subtree(tree, subset).internal_nodes() ]) for subset in subsets]
        for terrace_num, terrace_subtrees in enumerate(terrace_subtree_list):
            #if all subtrees matched, terrace is same
            if terrace_subtrees == this_tree_subtrees:
                out_stream.write('tree %d is on terrace %d\n' % (tree_num, terrace_num))
                terrace_size[terrace_num] += 1
                break
        else:
            #get here if went through all terraces with no match
            out_stream.write('tree %d is on a new terrace\n' % tree_num)
            terrace_size[len(terrace_subtree_list)] = 1
            terrace_subtree_list.append(this_tree_subtrees)

    messages.write(' done.')
    out_stream.write('%d terraces found' % len(terrace_size))


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


def open_tree_viewer(viewer_command, treefile, tree_object=None):
    '''If only a filename is passed in, try to open it.  If a tree object is passed, 
    write it to that filename and then try to open it.
    Note that the file extension of the treefile may matter to whatever is trying to open it
    '''
    if tree_object:
        with open(treefile, 'w') as outtree:
            outtree.write('%s\n' % tree_object)
    viewer_command =  shlex.split(viewer_command)
    viewer_command.append(treefile)
    try:
        #launch the command in the background
        proc = subprocess.Popen(viewer_command)
    except OSError:
        sys.stderr.write('Problem opening external tree viewer!\n')
        pass


########################################

parser = ArgumentParser(description='Perform various analyses related to phylogenetic terraces. Invoke script without arguments to start Tk GUI.')

in_group = parser.add_argument_group('Input Files')

in_group.add_argument('--alignment-file', default=None, help='nexus alignment including charsets to be used to determine character partition')

in_group.add_argument('--parent-tree-file', default=None, help='single parent tree to be analyzed (nexus or newick)')

in_group.add_argument('--treefiles-to-assign', nargs="*", default=None, help='trees to be assigned to one or more terraces (nexus or newick)')

in_group.add_argument('--subset-file', default=None, help='file with lines indicating sets of taxa represented in various partition subsets (created by --coverage preprocessing option)')

in_group.add_argument('--subtree-file', default=None, help='tree file containing the subtrees induced by the coverage subsets')

in_group.add_argument('--triplet-file', default=None, help='tab or space delimited triplet file, with ingroup ingroup outgroup (created by --triplets preprocessing option)')

preprocess = parser.add_argument_group('Preprosessing steps to perform on input files.  \nGeneral workflow would be --coverage, --display and --triplets, with each creating output consumed by following steps')

preprocess.add_argument('-c', '--coverage', action='store_true', default=False, help='compute the taxon coverage matrix (aka subsets file, requires --alignment-file)')

preprocess.add_argument('-d', '--display', action='store_true', default=False, help='print the subtrees displayed by the input tree with the input subsets (requires --subset-file and --tree-files)')

preprocess.add_argument('-t', '--triplets', action='store_true', default=False, help='Output arbitrary rooted taxon triples defining each edge in a set of treefiles (requires --tree-files)')

analyses = parser.add_argument_group('Analyses to be performed on files created by preprocessing')

#analyses.gui_options = {'start_hidden':True}

analyses.add_argument('-p', '--parents', action='store_true', default=False, help='compute the number of parent trees given a triplets file (requires --triplet-file')

analyses.add_argument('-b', '--build', action='store_true', default=False, help='compute the BUILD tree from a triplet file (requires --triplet-file)')

analyses.add_argument('-s', '--strict', action='store_true', default=False, help='compute a strict consensus tree from a triplet file (requires --triplet-file)')

analyses.add_argument('-l', '--list-terraces', action='store_true', default=False, help='take a set of trees and assign them to terraces (requires --subset-file and --tree-files)')

parser.add_argument('--simulate-coverage', type=float, nargs=3, default=None, help='simulate coverage matrices, using 3 values, #taxa #loci coverage')

parser.add_argument('--profile', action='store_true', default=False, help='profile the given functionality')

parser.add_argument('--open-tree-viewer', action='store_true', default=False, help='attempt to open an external tree viewer for build or strict consensus trees')

parser.add_argument('--verbose', action='store_true', default=False, help='spit out extra information for debugging purposes')

stdout_writer = MultiWriter(sys.stdout.write)
stderr_writer =  MultiWriter(sys.stderr.write)

#if no arguments are passed, try to start the tkinter gui
if len(sys.argv) == 1:
    try:
        from Tkinter import *
        from tkarg.tkinterutils import *
        #from ttk import *
    except ImportError:
        stderr_writer.write('%s\n' % parser.format_help())
        stderr_writer.write('\nUnable to import GUI componenets.  Use command line options.\n\n'.upper())
        sys.exit()
    
    preprocess.GUI_IGNORE = True

    tk_root = Tk()
    tk_gui = ArgparseGui(parser, tk_root, width=1200, height=720, destroy_when_done=False, progress_bar=False)

    #indicate dependencies between options to the gui, which will cause the dependent options to be disabled (greyed out) until the
    #dependency has been entered.  Options may have multiple dependencies.
    
    tk_gui.register_dependencies({'--triplet-file':['--build', '--parents', '--strict'], 
                                '--alignment-file':'--coverage', 
                                '--subset-file':['--display', '--list-terraces'],
                                '--parent-tree-file':'--display', 
                                '--treefiles-to-assign':'--list-terraces'})
    

    stderr_writer.add_streams(tk_gui.write_to_status)
   
    opt = tk_gui.get_option('--subset-file')
    opt.add_save_and_callback_button('COMPUTE', print_subsets, tk_gui.get_option('--alignment-file').file_count, tk_gui.get_option('--alignment-file').var, messages=stderr_writer)
    
    opt = tk_gui.get_option('--subtree-file')
    opt.add_save_and_callback_button('COMPUTE', print_displayed_subtrees, tk_gui.get_option('--subset-file').file_count, tk_gui.get_option('--parent-tree-file').var, 
            tk_gui.get_option('--subset-file').var, messages=stderr_writer)
    
    opt = tk_gui.get_option('--triplet-file')
    opt.add_save_and_callback_button('COMPUTE', calculate_triplets, tk_gui.get_option('--subtree-file').file_count, tk_gui.get_option('--subtree-file').var, messages=stderr_writer)
    
    #wait for the user to interact with the gui and press RUN button
    tk_root.mainloop()
   
    if tk_gui.cancelled:
        sys.exit('cancelled ...')
   
    #make the command line list generated by the gui that will be parsed by argparse exactly as if it were
    #typed on the command line
    options = parser.parse_args(tk_gui.make_commandline_list())
    prof = cProfile.Profile() if options.profile else None
    
    #set up a window to hold various result output
    results_window = ResultsWindow(tk_root, 300, 20)
    text = results_window.add_text_pane(0, 0)
    #to use the MultiWriter class the write function for the output window can only take one 
    #argument, hence the use of partial here
    from functools import partial
    stdout_writer.add_streams(partial(text.insert, END))

else:
    tk_root = None
    options = parser.parse_args()

    prof = cProfile.Profile() if options.profile else None
    #these pre-processing steps are done differently with the GUI
    if options.coverage:
        if not options.alignment_file:
            sys.exit('alignment file (-a) must be supplied to determine taxon coverage')
        profile_wrapper(print_subsets, prof, stdout_writer, options.alignment_file, messages=stderr_writer)

    if options.display:
        if not options.subset_file and options.parent_tree_file:
            sys.exit('must specify both subset file (-s) and tree file (--tree-files) to print displayed subtrees')
        profile_wrapper(print_displayed_subtrees, prof, stdout_writer, options.parent_tree_file, options.subset_file, messages=stderr_writer)

    if options.triplets:
        if not options.subtree_file:
            sys.exit('subtree file (--subtree-file) must be specified to calculate triplets')
        profile_wrapper(calculate_triplets, prof, stdout_writer, options.subtree_file, stderr_writer)

#for debugging, so gui doesn't actually need to be manipulated
if tk_root:
    #if not options.subset_file:
    #    options.subset_file = 'subsets'
    if not options.triplet_file:
        options.triplet_file = 'triplets'
    #if not options.tree_files:
    #    options.tree_files = ['subtrees.tre']
    #options.build = True
    #options.strict = True
    pass

if options.simulate_coverage:
    taxa = int(options.simulate_coverage[0])
    loci = int(options.simulate_coverage[1])
    cov = options.simulate_coverage[2]

    sim_tree = uniform_pure_birth(TaxonSet(['t%d' % num for num in xrange(taxa)]))
    sim_tree.write_to_path('parent.tre', schema='nexus')

    mat = CoverageMatrix()
    mat.fill_random(taxa, loci, cov, lambda x, y, c: c, reference_taxon=True)
    mat.print_subset_vectors()
    #exponential pdf is lam * e^(lam*x)
    #where scale param is 1/lam=mean
    #exp_scale = 1.0
    from functools import partial
    from numpy.random import exponential
    #mat.fill_random_locus_func(taxa, loci, partial(exponential, cov), min_coverage=0.2, max_coverage=0.8, reference_taxon=True)
    #mat.fill_random_taxon_func(taxa, loci, partial(exponential, cov), min_coverage=0.2, max_coverage=0.8, reference_taxon=True)
    #mat.print_subset_vectors()
    #profile_wrapper(mat.test_decisiveness, prof)

    prob = 1.0 - (taxa - 2.0) * (1.0 - cov**3)**loci
    #print 'Min prob decisive for tree: %f' % prob
    #stderr_writer.write('Min prob decisive for tree: %f\n' % prob)

if options.subset_file:
    
    subsets = read_subset_file(options.subset_file)
    have_subsets = True
    
    mat = CoverageMatrix()
    mat.fill_from_subsets(subsets)

    loci_per_taxon =  mat.loci_per_taxon()
    mat.calculate_statistics()
    
    if tk_root:
        bar_canvas = results_window.add_canvas_pane(0, 1)
        mat.draw_barplot(bar_canvas, loci_per_taxon, 0, 0, bar_canvas.winfo_reqwidth(), bar_canvas.winfo_reqheight())
        
        matrix_canvas = results_window.add_canvas_pane(1, 0, column_span=2)
        mat.draw_matrix_graphic(matrix_canvas, mat.taxa, 0, 0, matrix_canvas.winfo_reqwidth(), matrix_canvas.winfo_reqheight())
        
        def sort_and_redraw():
            taxa = list(mat.taxa)
            def cells(tax):
                return sum(mat.per_taxon_presence_absence[tax])
            taxa.sort(key=cells, reverse=True)
            mat.draw_matrix_graphic(matrix_canvas, taxa, 0, 0, matrix_canvas.winfo_reqwidth(), matrix_canvas.winfo_reqheight())

        button = Button(matrix_canvas, text='SORT', command=sort_and_redraw)
        matrix_canvas.create_window((0, 0), window=button, height=button.winfo_reqheight(), width=button.winfo_reqwidth(), anchor='nw')
        
    else:
        out_trans = ['-', 'X']
        for tax, cov in sorted(mat.per_taxon_presence_absence.items()):
            sys.stderr.write('%30s\t%s\n' % (tax, ''.join([out_trans[c] for c in cov])))

while True:
    if options.build:
        if not options.triplet_file:
            sys.exit('triplet file (-t) must be supplied to make BUILD tree')
        
        build_tree = profile_wrapper(make_build_tree, prof, stdout_writer, options.triplet_file, messages=stderr_writer, verbose=options.verbose)
       
        if options.open_tree_viewer:
            open_tree_viewer(tree_viewer_command, 'build.tre', build_tree)
        
    if options.strict:
        if not options.triplet_file:
            sys.exit('triplet file (-t) must be supplied to make strict consensus tree')
        
        use_threads = False
        if use_threads and tk_gui:
            #threading is sort of working here, but needs a lot more work
            #tk_gui.progress_bar.start()

            class ThreadedTask(threading.Thread):
                def __init__(self, queue, func, *args, **kwargs):
                    threading.Thread.__init__(self)
                    self.queue = queue
                    self.func = func
                    self.args = args
                    self.kwargs = kwargs
                def run(self):
                    result = self.func(*self.args, **self.kwargs)
                    self.queue.put(result)

            tk_gui.queue_thread(ThreadedTask(tk_gui.queue, make_strict_tree, stdout_writer, options.triplet_file, messages=stderr_writer, verbose=options.verbose))
            tk_root.mainloop()
            tk_root.mainloop()
            strict_tree = tk_gui.result
            
            #tk_gui.progress_bar.stop()
        else: 
            strict_tree = profile_wrapper(make_strict_tree, prof, stdout_writer, options.triplet_file, messages=stderr_writer, verbose=options.verbose)
        
        if options.open_tree_viewer:
            open_tree_viewer(tree_viewer_command, 'strict.tre', strict_tree)
        
    if options.parents:
        if not options.triplet_file:
            sys.exit('triplet file (-t) must be supplied to count parent trees')
        profile_wrapper(count_trees_on_terrace, prof, stdout_writer, options.triplet_file, messages=stderr_writer)
        if tk_root:
            tk_root.mainloop()

    if options.list_terraces:
        if not options.subset_file or not options.treefiles_to_assign:
            sys.exit('must specify both subset file (-s) and --trees-to-assign to assign trees to terraces')
        profile_wrapper(assign_to_terraces, prof, stdout_writer, options.treefiles_to_assign, options.subset_file, messages=stderr_writer)

    if tk_root:
        tk_root.mainloop()
        #update the options namespace to match the current state of the gui
        options = parser.parse_args(tk_gui.make_commandline_list())
    else:
        break

if prof:
    output_profile(prof)


