Scripts for analyzing phylogenetic terraces.
Derrick J. Zwickl

INSTALLATION:

Obtain the source and example files in one of the following ways:
checkout source:
git clone https://github.com/zwickl/terraphy.git
svn checkout https://github.com/zwickl/terraphy/trunk terraphy
download:
https://github.com/zwickl/terraphy/archive/master.zip

enter the terraphy directory, then the typical python install
procedure of 
python setup.py install 
should work (you might need to add sudo to the start of that line).  

Dependencies should automatically be installed with the above command.

If you get errors with that command, you might need to update your python 
"setuptools" package, i.e.:
easy_install setuptools
(again, you might need to add sudo to the start of that line)

Or you can just run things in the example directories, as demonstrated in
the test.sh sample scripts.  (although you'll need to install the dendropy
dependency yourself)

If you install it globally you may need to put the script install location 
into your PATH environment variable, or copy the script files (in the 
scripts directory) somewhere else in your path.

Dependencies:
-The excellent Dendropy library:

Sukumaran J, Holder MT. 2010. DendroPy: a Python library for phylogenetic 
    computing. Bioinformatics. 26(12):1569-1571.
http://pythonhosted.org/DendroPy/


