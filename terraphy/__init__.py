import sys
import os

###############################################################################
## Populate the 'terraphy' namespace

from terraphy import triplets

###############################################################################
## PACKAGE METADATA

__project__ = "terraphy"
__version__ = "1.0"

try:
    try:
        __homedir__ = __path__[0]
    except AttributeError:
        __homedir__ = os.path.dirname(os.path.abspath(__file__))
    except IndexError:
        __homedir__ = os.path.dirname(os.path.abspath(__file__))
except OSError:
    __homedir__ = None
except:
    __homedir__ = None

#__revision__ = vcsinfo.Revision(repo_path=__homedir__)
__author__ = "Derrick Zwickl"
__copyright__ = "Copyright 2014 Derrick Zwickl"
__license__ = """
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program. If not, see <http://www.gnu.org/licenses/>.
"""
PACKAGE_VERSION = __version__ # for backwards compatibility (with sate)

def description():
    if __revision__.is_available:
        revision_text = " (%s)" % str(__revision__)
    else:
        revision_text = ""
    return "%s %s%s" % (__project__, __version__, revision_text)

if __name__ == "__main__":
    sys.stdout.write("%s\n" % description())


