#!/usr/bin/env python
import sys
import re
from argparse import ArgumentParser

sys.path.append("../")

from terraphy.triplets import find_triplets_defining_edges_descending_from_node

import dendropy

#for dendropy 4 compatability
try:
    from dendropy.error import DataError
except ImportError:
    from dendropy.utility.error import DataError

parser = ArgumentParser(description='Read newick or nexus input trees from file or stdin and output arbitrary taxon triples defining each edge')

parser.add_argument('treefiles', nargs='*', default=[], help='nexus or newick treefile(s) to read (omit for stdin)')

options = parser.parse_args()


intrees = dendropy.TreeList()
if not options.treefiles:
    sys.stderr.write('NOTE: reading trees from stdin\n')
    trees = sys.stdin.read()
    #try two input formats
    try:
        intrees.extend(dendropy.TreeList.get_from_string(trees, "nexus"))
    except DataError:
        intrees.extend(dendropy.TreeList.get_from_string(trees, "newick"))
else:
    for tf in options.treefiles:
        #try two input formats
        try:
            intrees.extend(dendropy.TreeList.get_from_path(tf, "nexus"))
        except DataError:
            intrees.extend(dendropy.TreeList.get_from_path(tf, "newick"))
        except ValueError:
            sys.stderr.write('NOTE: ValueError reading from file %s, ' % tf)
        except AttributeError:
            sys.stderr.write('NOTE: AttributeError reading from file %s, ' % tf)

sys.stderr.write('read %d trees\n' % len(intrees))

triplets = []
for tree in intrees:
    treestr = tree.as_newick_string()
    treestr = re.sub('([A-Za-z_.-]+)', '"\\1"', treestr)
    triplets.extend(find_triplets_defining_edges_descending_from_node(eval(treestr)))

#collect all labels included in any triplet, which should be the actual number 
#included in the input trees
all_taxa = set()
for trip in triplets:
    for tax in trip:
        all_taxa.add(tax)

print ' '.join([tax for tax in sorted(all_taxa)])
print '\n'.join([' '.join(trip) for trip in triplets])

