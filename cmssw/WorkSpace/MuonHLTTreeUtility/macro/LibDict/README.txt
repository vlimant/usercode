This collection of files generates the libraries and dictionaries needed by root to handle data types of the forms
map <int, vector<int> >
map <int, vector<double> >

Execution:
From this directory, execute the Makefile with the make command.  The output is the libdict file libTreeLib.so.  For executing the macros that make use of the data types listed above, execute either in the rootlogon.C file or from the command line

gSystem->Load("libTreeLib.so");

before loading the actual macro.