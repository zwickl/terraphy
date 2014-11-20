#!/usr/bin/env python
import sys
import re
from argparse import ArgumentParser

sys.path.append("../")

from terraphy.triplets import find_triplets_defining_edges_descending_from_node

try:
    import dendropy
except ImportError:
    sys.exit('problem importing dendropy package - it is required')

from terraphy.dendroutils import dendropy_read_treefile

parser = ArgumentParser(description='Read newick or nexus input trees from file or stdin and output arbitrary taxon triples defining each edge')

parser.add_argument('treefiles', nargs='*', default=[], help='nexus or newick treefile(s) to read (omit for stdin)')

options = parser.parse_args()

intrees = dendropy_read_treefile(options.treefiles, suppress_edge_lengths=True)
sys.stderr.write('read %d trees\n' % len(intrees))

triplets = []
for tnum, tree in enumerate(intrees):
    if hasattr(tree, 'as_newick_string'):
        treestr = tree.as_newick_string()
    else:
        treestr = tree.as_string(schema='newick', suppress_edge_lengths=True, suppress_rooting=True)
    #eliminate any singletons
    treestr = re.sub('([A-Za-z_.-]+)', '"\\1"', treestr)
    treestr = re.sub(';', '', treestr)
    #evaluate newick string as set of nested tuples
    triplets.extend(find_triplets_defining_edges_descending_from_node(eval(treestr)))
    sys.stderr.write('after tree %d: %d triplets\n' % (tnum, len(triplets)))

#collect all labels included in any triplet, which should be the actual number 
#included in the input trees
all_taxa = set()
for trip in triplets:
    for tax in trip:
        all_taxa.add(tax)

sys.stdout.write('%s\n' % ' '.join([tax for tax in sorted(all_taxa)]))
sys.stdout.write('%s\n' % '\n'.join([' '.join(trip) for trip in triplets]))

