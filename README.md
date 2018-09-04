# Generators
For playing w/ PYTHIA8, Herwig++, etc.

From top dir, do "source setGeneratorsPath.sh". Make will not work w/o this

cd int PYTHIA8
run make
this requires a few environmental variables
 * check that "root-config --help" returns something meaningful
 * check that "pythia8-config --help" returns something meaningful
 * check that "lhapdf-config --help" returns something meaningful
 * check that "$FASTJETPATH/bin/fastjet-config --help" returns something meaningful
 
If "root-config --help" does not return meaningful output, you really have to install/properly configure root to proceed

If "pythia8-config --help" or "lhapdf-config --help" or "$FASTJETPATH/bin/fastjet-config --help" does not return meaningful output, you can proceed by
 * Replacing the line "all: ..." in Makefile removing "bin/pythia8CUETP8M2T4.exe"

Then try running make again
if make is successful, try the statistical comparison example script:

./bin/statisticalComparison.exe

this will tell you the usage is:

Usage: ./bin/statisticalComparison.exe <flatPthatFileName> <stagPthatFileName>

try again with correct inputs. Use this macro as an analysis template