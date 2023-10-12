#!/bin/bash

#module load GCC/8.3.0 OpenMPI/3.1.4 ROOT/6.20.04-Python-3.7.4
#source /opt/ebsofts/ROOT/6.20.04-foss-2019b-Python-3.7.4/bin/thisroot.sh

#Analyse list of root files 
rootFilesList="./rootFile.list"
outHistF="./hist.root"

#Or analyse single root file 
inRootFiles="../scratch/cosmique_gamma_hadron_generator/out/00010/cosmique_gamma_hadron_generator.root"
outHistSingleF="./hist_one.root"

make -f Makefilecpv clean; make -f Makefilecpv runcpv;

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d  : single root file"
    echo " [0] -l  : list of root files"
    echo " [0] -h  : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	./runcpv 1 $inRootFiles $outHistSingleF
    elif [ "$1" = "-l" ]; then
	./runcpv 0 $rootFilesList $outHistF
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak "I have done"
