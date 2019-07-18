To install:


in RNAstructure dir:

make python_setup

get an error 
change newly generated .so file name to _RNAstructure_wrap.so

make python_setup again.

the files in exe folder is the python modules.

move exe folder and data_tables folder to a RNAFolder.

for proper import to work:

add following to .bash_profile


PYTHONPATH+="/Library/Frameworks/Python.framework/Versions/3.7/lib/python3.7/site-packages/RNAstructure/exe"
export PYTHONPATH
export DATAPATH=/Library/Frameworks/Python.framework/Versions/3.7/lib/python3.7/site-packages/RNAstructure/data_tables/
