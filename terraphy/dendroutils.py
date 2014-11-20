#/usr/bin/env python
import sys

#DENDROPY PACKAGE
try:
    from dendropy import TreeList
    
    #for dendropy 4 compatability
    try:
        from dendropy.error import DataError
    except:
        from dendropy.utility.error import DataError

except ImportError:
    sys.exit('problem importing dendropy package - it is required')


def dendropy_read_treefile(treefiles, quiet=False, **kwargs):
    intrees = TreeList()
    if not treefiles:
        if not quiet:
            sys.stderr.write('NOTE: reading trees from stdin\n')
        trees = sys.stdin.read()
        #try two input formats
        try:
            intrees.extend(TreeList.get_from_string(trees, "nexus", case_sensitive_taxon_labels=True, preserve_underscores=True, **kwargs))
        except DataError:
            intrees.extend(TreeList.get_from_string(trees, "newick", case_sensitive_taxon_labels=True, preserve_underscores=True, **kwargs))
    else:
        for tf in treefiles:
            #try two input formats
            try:
                if not quiet:
                    sys.stderr.write('Reading file %s in newick format ...\n' % tf)
                intrees.extend(TreeList.get_from_path(tf, "newick", case_sensitive_taxon_labels=True, preserve_underscores=True, **kwargs))
            except DataError:
                if not quiet:
                    sys.stderr.write('Reading file %s in nexus format ...\n' % tf)
                intrees.extend(TreeList.get_from_path(tf, "nexus", case_sensitive_taxon_labels=True, preserve_underscores=True, **kwargs))
            except ValueError:
                sys.stderr.write('NOTE: ValueError reading from file %s, ' % tf)
            except AttributeError:
                sys.stderr.write('NOTE: AttributeError reading from file %s, ' % tf)

    return intrees


