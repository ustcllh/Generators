#!/bin/bash
filepath="/eos/cms/store/group/phys_heavyions/cmcginn/PYTHIA8_ak4GenSkim_20180924/"
./bin/DijetImbalanceRatio_pythia.exe $(find $filepath -name *Pythia8*All*)
mv output/statComp*.root output/CUETP8M1_CMS.root
