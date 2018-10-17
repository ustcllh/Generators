#!/bin/bash
#filepath="/eos/cms/store/group/phys_heavyions/cmcginn/PYTHIA6_ak4GenSkim_20180916/"
#use smearing data
filepath="/eos/cms/store/group/phys_heavyions/cmcginn/PYTHIA6_ak4GenSkim_20180924/"
./bin/DijetImbalanceRatio_pythia.exe $(find $filepath -name *Pythia6*All*)
mv output/statComp*.root output/PYTHIA6_CMS.root
