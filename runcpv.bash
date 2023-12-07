#!/bin/bash

#module load GCC/8.3.0 OpenMPI/3.1.4 ROOT/6.20.04-Python-3.7.4
#source /opt/ebsofts/ROOT/6.20.04-foss-2019b-Python-3.7.4/bin/thisroot.sh

#Analyse list of root files 
rootFilesList="./rootFile120km.list"
outHistF="./hist120km.root"

#Or analyse single root file 
inRootFiles="./cosmique_gamma_hadron_generator_corsika.root"
outHistSingleF="./hist_corsika.root"

make -f Makefilecpv clean; make -f Makefilecpv runcpv;

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d                          : single root file"
    echo " [0] -l                          : list of root files"
    echo " [0] -ls                         : list (short) of root files"
    echo " [0] -fl                         : calculate fluxes"
    echo " [0] -proton_diff                : proton diff."
    echo " [0] -gamma_diff                 : gamma diff."
    echo " [0] -gamma                      : gamma on axis"
    echo " [0] -ele_pos                    : electron"
    echo " [0] -proton_diff_cghg_vs_simtel : proton diff. cghg vs simtel"
    echo " [0] -gamma_diff_cghg_vs_simtel  : gamma diff. cghg vs simtel"
    echo " [0] -gamma_cghg_vs_simtel       : gamma on axis cghg vs simtel"
    echo " [0] -ele_pos_cghg_vs_simtel     : electron cghg vs simtel"
    echo " [0] -h                          : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	./runcpv 1 $inRootFiles $outHistSingleF
    elif [ "$1" = "-l" ]; then
	time ./runcpv 0 $rootFilesList $outHistF | tee log
    elif [ "$1" = "-ls" ]; then
	time ./runcpv 0 "./rootFile120km_short.list" $outHistF | tee log
    elif [ "$1" = "-fl" ]; then
	time ./runcpv 2 | tee log
    elif [ "$1" = "-proton_diff" ]; then
	rootFilesList="./rootFile120km.list"
	outHistF="./hist120km_proton_diff.root"
	time ./runcpv 3 $rootFilesList $outHistF proton_diff 100000000
    elif [ "$1" = "-gamma_diff" ]; then
	rootFilesList="./rootFile120km.list"
	outHistF="./hist120km_gamma_diff.root"
	time ./runcpv 3 $rootFilesList $outHistF gamma_diff 100000000
    elif [ "$1" = "-gamma" ]; then
	rootFile="../scratch/cosmique_gamma_hadron_generator/120km_onaxis/out/00000/cosmique_gamma_hadron_generator_corsika_onaxis.root"
	outHistF="./hist120km_gamma.root"
	time ./runcpv 4 $rootFile $outHistF gamma 100000000
    elif [ "$1" = "-ele_pos" ]; then
	rootFilesList="./rootFile120km.list"
	outHistF="./hist120km_ele_diff.root"
	time ./runcpv 3 $rootFilesList $outHistF ele_pos 100000000
    elif [ "$1" = "-proton_diff_cghg_vs_simtel" ]; then
	rootFilesList="./rootFile120km.list"
	outHistF="./hist120km_proton_diff_cghg_vs_simtel.root"
	time ./runcpv 5 $rootFilesList $outHistF proton_diff 100000000
    elif [ "$1" = "-gamma_diff_cghg_vs_simtel" ]; then
	rootFilesList="./rootFile120km.list"
	outHistF="./hist120km_gamma_diff_cghg_vs_simtel.root"
	time ./runcpv 5 $rootFilesList $outHistF gamma_diff 10000000
    elif [ "$1" = "-gamma_cghg_vs_simtel" ]; then
	rootFile="../scratch/cosmique_gamma_hadron_generator/120km_onaxis/out/00000/cosmique_gamma_hadron_generator_corsika_onaxis.root"
	outHistF="./hist120km_gamma_cghg_vs_simtel.root"
	time ./runcpv 6 $rootFile $outHistF gamma 10000000
    elif [ "$1" = "-ele_pos_cghg_vs_simtel" ]; then
	rootFilesList="./rootFile120km.list"
	outHistF="./hist120km_ele_diff_cghg_vs_simtel.root"
	time ./runcpv 5 $rootFilesList $outHistF ele_pos 10000000	
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak "I have done"
