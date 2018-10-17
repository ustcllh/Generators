#!/bin/bash
filepath="/eos/cms/store/group/phys_heavyions/cmcginn/PYTHIA8Tunes_20180926/"
./bin/DijetImbalanceRatio.exe $(find $filepath -name *Flat*CP5*) $(find $filepath -name *Pthat*CP5*)
mv output/statComp*.root output/CP5.root
./bin/DijetImbalanceRatio.exe $(find $filepath -name *Flat*CUETP8M1*) $(find $filepath -name *Pthat*CUETP8M1*)
mv output/statComp*.root output/CUETP8M1.root
./bin/DijetImbalanceRatio.exe $(find $filepath -name *Flat*CUETP8M2T4*) $(find $filepath -name *Pthat*CUETP8M2T4*)
mv output/statComp*.root output/CUETP8M2T4.root
./bin/DijetImbalanceRatio.exe $(find $filepath -name *Flat*Vanilla*) $(find $filepath -name *Pthat*Vanilla*)
mv output/statComp*.root output/Vanilla.root
