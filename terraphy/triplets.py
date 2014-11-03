#!/usr/bin/env python
import sys
import re
from argparse import ArgumentParser
from collections import Iterable


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
    '''

    if isinstance(components, str):
        #shouldn't be getting here because the function won't be recursively called with a leaf
        sys.exit('find_triplets_defining_edges_descending_from_node() called with leaf?')
    if len(components) == 1:
        return []

    toReturn = []
    for in_subtree in components:
        if len(in_subtree) > 1 and not isinstance(in_subtree, str):
            for out_subtree in components:
                if out_subtree != in_subtree:
                    gen = flattened_array_generator(in_subtree[0], sys.maxint)
                    ingroup1 = next(gen)
                    gen = flattened_array_generator(in_subtree[1], sys.maxint)
                    ingroup2 = next(gen)
                    gen = flattened_array_generator(out_subtree, sys.maxint)
                    outgroup = next(gen)
                    #ingroup1 = [f for f in flattened_array_generator(in_subtree[0], sys.maxint)][0]
                    #ingroup2 = [f for f in flattened_array_generator(in_subtree[1], sys.maxint)][0]
                    #outgroup = [f for f in flattened_array_generator(out_subtree, sys.maxint)][0]
                    toReturn.append([ingroup1, ingroup2, outgroup])
                    #break so that different triplets are not chosen that only differ by
                    #outgroup, in the case of polytomies
                    break

    for comp in components:
        if len(comp) > 1 and not isinstance(comp, str):
            toReturn.extend(find_triplets_defining_edges_descending_from_node(comp))

    return toReturn


if __name__ == "__main__":
    import doctest
    doctest.testmod()

