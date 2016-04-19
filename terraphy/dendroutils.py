#/usr/bin/env python
import sys
import os

#DENDROPY PACKAGE
try:
    from dendropy import TreeList
    
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


def dendropy_read_treefile(treefiles, quiet=False, **kwargs):
    out_stream = kwargs.pop('writer', sys.stderr)
    intrees = TreeList()
    if not treefiles:
        if not quiet:
            sys.stderr.write('NOTE: reading trees from stdin\n')
        trees = sys.stdin.read()
        #try two input formats
        try:
            intrees.extend(TreeList.get_from_string(trees, "nexus", case_sensitive_taxon_labels=True, preserve_underscores=True, **kwargs))
        except (DataParseError, NexusReader.NotNexusFileError) as e:
            sys.stderr.write('%s\n' % e.message)
            intrees.extend(TreeList.get_from_string(trees, "newick", case_sensitive_taxon_labels=True, preserve_underscores=True, **kwargs))
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
                intrees.extend(TreeList.get_from_path(tf, "nexus", case_sensitive_taxon_labels=True, preserve_underscores=True, **kwargs))

            #except (DataParseError, dendropy.dataio.nexusreader.NotNexusFileError) as e:
            except (DataParseError, NexusReader.NotNexusFileError, AttributeError) as e:
                try:
                    if not quiet:
                        out_stream.write('Reading file %s in newick format ...\n' % tf)
                    intrees.extend(TreeList.get_from_path(tf, "newick", case_sensitive_taxon_labels=True, preserve_underscores=True, **kwargs))
                except (DataParseError, Tokenizer.UnexpectedEndOfStreamError, AttributeError)  as e:
                    if not quiet:
                        sys.stderr.write('%s\n' % e.message)
                        sys.exit('Could not read file %s in nexus or newick  format ...\n' % tf)
    return intrees


