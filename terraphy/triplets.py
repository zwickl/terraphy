#!/usr/bin/env python
import sys
import re
from argparse import ArgumentParser
from collections import Iterable
from itertools import combinations


def flattened_array_generator(array, level=1, reverse=False):
    '''Generator to be used in flatten_array function, or by itself.
    Difference is that this doesn't create the flattended list, as
    flatten_array does, it just yields elements as if it had.
    Favor this in iteration, obviously.
    '''
    if isinstance(array, Iterable) and not isinstance(array, str):
        if reverse:
            array = array[::-1]

        for toplvl in array:
            if level == 0:
                yield toplvl
            else:
                for sub in flattened_array_generator(toplvl, level=level - 1, reverse=reverse):
                    yield sub
    else:
        yield array

def ingroup_pair(ingroup, all_pairs=False):
    '''generator that either yields all pairs of elements of in_subtree,
    or just the first element paired with each of the others
    for use in find_triplets_defining_edges_descending_from_node
    '''

    if all_pairs:
        for ingroup1, ingroup2 in combinations(in_subtree, 2):
            yield ingroup1, ingroup2
    else:
        for ingroup2 in ingroup[1:]:
            yield ingroup[0], ingroup2


def find_triplets_defining_edges_descending_from_node(components):
    '''Recursively find triplets that "define" each of the edges descending from a node.  
    At node Z, edges to define are Z-X and Z-Y.  A, B, C and D may be leaves or subtrees:

    A   B C   D   
     \ /   \ /
      X     Y
       \   /
         Z

    To define Z-X, we need two "ingroup" taxa, each arbitrarily drawn from different descendents 
    of X, in this case from A and B.  If X had > 2 descendents, pick any two.  We also need an 
    outgroup taxon, which for Z-X would be any taxon in clade Y.  Likewise, to define Z-Y, need 
    one taxon drawn from each of clades C, D and X.
    
    Input - components is just a set of nested tuples (or lists), equivalent to the newick 
    representation of clade Z.
    Output - list of ingroup ingroup outgroup triplets, one for each edge defined.  Specific taxa 
    appearing in triplets are arbitray for many edges, although are consistent based on the 
    rotation of clades in a given newick string. 

    With polytomies things are more complicated
    A  B  E  
     \ | / C  D
       X   | /   
       |   Y     
        \/
        Z
    Given the above, to define X-Z we need triplets
    AB|C
    AE|C
    i.e. it isn't enough to arbitrarily pick two of A, B, E
    nor do all combinations need to be generated
    generator ingroup_pair is used to get appropriate ingroup elements
    '''

    if isinstance(components, str) or len(components) == 1:
        return []

    if len(components) > 2:
        sys.stderr.write("WARNING: making triplets from a tree with polytomies!\n")
        poly = True
    else:
        poly = False

    toReturn = []
    for in_subtree in components:
        if len(in_subtree) > 1 and not isinstance(in_subtree, str):
            for out_subtree in components:
                if out_subtree != in_subtree:
                    trip = []

                    all_pairs = False

                    #for in_subtree1, in_subtree2 in combinations(in_subtree, 2):
                    for in_subtree1, in_subtree2 in ingroup_pair(in_subtree):
                        ingroup1 = sorted([ tax for tax in flattened_array_generator(in_subtree1, sys.maxsize) ])[0]
                        ingroup2 = sorted([ tax for tax in flattened_array_generator(in_subtree2, sys.maxsize) ])[0]
                        #print in_subtree, in_subtree1, in_subtree2, ingroup1, ingroup2, out_subtree
                    
                        alpha_names = True
                        if alpha_names:
                            #Choose the alphabetically first taxon from each subtree - this will take a bit of extra
                            #time to sort, but will allow the later identification and removal of triplets from different 
                            #trees that define the same edge but would otherwise appear different due to different 
                            #arbitrarily chosen taxa


                            #ingroup1 = sorted([ tax for tax in flattened_array_generator(in_subtree[0], sys.maxsize) ])[0]
                            #ingroup2 = sorted([ tax for tax in flattened_array_generator(in_subtree[1], sys.maxsize) ])[0]
                            
                            outgroup = sorted([ tax for tax in flattened_array_generator(out_subtree, sys.maxsize) ])[0]
                            trip = [ingroup1, ingroup2, outgroup]
                            #the ordering of the ingroup taxa is also arbitrary, so make it alphabetical
                            #can't do this with polytomies though
                            if ingroup1 > ingroup2 and not poly:
                                trip[0], trip[1] = trip[1], trip[0]

                            #print trip
                        else:
                            for subtree in [in_subtree[0], in_subtree[1], out_subtree]:
                                gen = flattened_array_generator(subtree, sys.maxsize)
                                trip.append(next(gen))
                        

                        toReturn.append(trip)
                        #break so that different triplets are not chosen that only differ by
                        #outgroup, in the case of polytomies
                        #NO, don't think that's right.  If for example there were a node with two leaf descendnents and one
                        #clade, both leaves wouldn't be included in any triplet unless this loops multiple times and each
                        #is included in a triplet as an outgroup.
                        #break

    for comp in components:
        if len(comp) > 1 and not isinstance(comp, str):
            toReturn.extend(find_triplets_defining_edges_descending_from_node(comp))

    return toReturn


if __name__ == "__main__":
    import doctest
    doctest.testmod()

