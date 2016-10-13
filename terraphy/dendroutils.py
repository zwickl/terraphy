#/usr/bin/env python
import sys
import os

#DENDROPY PACKAGE
try:
    from dendropy import TreeList, Tree, Node
    
    #for dendropy 4 compatability
    try:
        from dendropy.error import DataError as DataParseError
    except:
        from dendropy.utility.error import DataParseError

    #import dendropy.dataio.nexusreader.NotNexusFileError as NotNexusFileError
    from dendropy.dataio.nexusreader  import NexusReader
    from dendropy.dataio.tokenizer  import Tokenizer
    #from  dendropy.dataio.nexusreader import NotNexusFileError
    #from  dendropy.dataio.nexusreader import NexusReaderError

    #this deals with changes in DendroPy 4
    try:
        from dendropy.calculate import treesplit
    except ImportError:
        from dendropy import treesplit

except ImportError as e:
    sys.stderr.write('%s\n' % e.message)
    sys.exit('problem importing dendropy package - it is required')


def compat_encode_bipartitions(tree, **kwargs):
    '''Convenience function dealing with different ways of encoding splits in DP4'''
    if hasattr(tree, "encode_bipartitions"):
        if 'delete_outdegree_one' in kwargs:
            val = kwargs.pop('delete_outdegree_one')
            kwargs['collapse_unrooted_basal_bifurcation'] = val
        tree.encode_bipartitions(**kwargs)
    elif not hasattr(tree, "bipartition_encoding") or not tree.bipartition_encoding:
        treesplit.encode_splits(tree, **kwargs)


def displayed_subtree(tree, labels):
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
            newtree.retain_taxa_with_labels(labels)
    else:
            newtree.retain_taxa(labels)

    #compat_encode_bipartitions now maps delete_outdegree_one to collapse_unrooted_basal_bifurcation in DP 4
    compat_encode_bipartitions(newtree, delete_outdegree_one=False)
    return newtree

def same_tree(reference_tree, test_tree):
    '''This is adapted from false_positives_and_negatives() in dendropy treecalc module,
    and just bails when it finds the first different split.
    Should handle polytomies fine.
    '''

    #ref_tset = compat_get_taxon_set(reference_tree)
    #test_tset = compat_get_taxon_set(test_tree)
    ref_tset = reference_tree.taxon_namespace
    test_tset =test_tree.taxon_namespace
    if ref_tset is not test_tset:
        raise TypeError("Trees have different TaxonSet objects: %s vs. %s" \
                % (hex(id(ref_tset)), hex(id(test_tset))))

    compat_encode_bipartitions(reference_tree)
    compat_encode_bipartitions(test_tree)

    '''
    #seems like this set comparison should be faster, but not really
    if isinstance(reference_tree.split_edges, dict):
        if set(reference_tree.split_edges.keys()) != set(test_tree.split_edges.keys()):
            return False
    elif set(reference_tree.split_edges) != set(test_tree.split_edges):
        return False
    '''
    for split in reference_tree.bipartition_encoding:
        if split not in test_tree.bipartition_encoding:
            return False

    for split in test_tree.bipartition_encoding:
        if split not in reference_tree.bipartition_encoding:
            return False

    return True




def dendropy_read_treefile(treefiles, quiet=False, preserve_underscores=False, **kwargs):
    out_stream = kwargs.pop('writer', sys.stderr)
    intrees = TreeList()
    if not treefiles:
        if not quiet:
            sys.stderr.write('NOTE: reading trees from stdin\n')
        trees = sys.stdin.read()
        #try two input formats
        try:
            intrees.extend(TreeList.get_from_string(trees, "nexus", case_sensitive_taxon_labels=True, preserve_underscores=preserve_underscores, **kwargs))
        except (DataParseError, NexusReader.NotNexusFileError) as e:
            sys.stderr.write('%s\n' % e.message)
            intrees.extend(TreeList.get_from_string(trees, "newick", case_sensitive_taxon_labels=True, preserve_underscores=preserve_underscores, **kwargs))
        except (DataParseError, Tokenizer.UnexpectedEndOfStreamError, AttributeError)  as e:
            if not quiet:
                sys.stderr.write('%s\n' % e.message)
                sys.exit('Could not read file %s in nexus or newick  format ...\n' % tf)
    else:
        for tf in treefiles:
            if not os.path.isfile(tf):
                out_stream.write('TreeFile %s  does not exist' % tf)
                sys.exit()

            #try two input formats
            try:
                if not quiet:
                    out_stream.write('Reading file %s in nexus format ...\n' % tf)
                intrees.extend(TreeList.get_from_path(tf, "nexus", case_sensitive_taxon_labels=True, preserve_underscores=preserve_underscores, **kwargs))

            #except (DataParseError, dendropy.dataio.nexusreader.NotNexusFileError) as e:
            except (DataParseError, NexusReader.NotNexusFileError, AttributeError) as e:
                try:
                    if not quiet:
                        out_stream.write('Reading file %s in newick format ...\n' % tf)
                    intrees.extend(TreeList.get_from_path(tf, "newick", case_sensitive_taxon_labels=True, preserve_underscores=preserve_underscores, **kwargs))
                except (DataParseError, Tokenizer.UnexpectedEndOfStreamError, AttributeError)  as e:
                    if not quiet:
                        sys.stderr.write('%s\n' % e.message)
                        sys.exit('Could not read file %s in nexus or newick  format ...\n' % tf)
    return intrees

