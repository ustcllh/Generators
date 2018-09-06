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


### set library path in setGeneratorsPath.sh
### use CMSSW_10_2_0
export PYTHIA8PATH="/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/pythia8/230-gnimlf4"
export LHAPDF6PATH="/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/lhapdf/6.2.1"
export FASTJETPATH="/cvmfs/cms.cern.ch/slc6_amd64_gcc700/external/fastjet/3.3.0"
 
If "root-config --help" does not return meaningful output, you really have to install/properly configure root to proceed

If "pythia8-config --help" or "lhapdf-config --help" or "$FASTJETPATH/bin/fastjet-config --help" does not return meaningful output, you can proceed by
 * Replacing the line "all: ..." in Makefile removing "bin/pythia8CUETP8M2T4.exe"

Then try running make again
if make is successful, try the statistical comparison example script:

./bin/statisticalComparison.exe

this will tell you the usage is:

Usage: ./bin/statisticalComparison.exe <flatPthatFileName> <stagPthatFileName>

try again with correct inputs. Use this macro as an analysis template
