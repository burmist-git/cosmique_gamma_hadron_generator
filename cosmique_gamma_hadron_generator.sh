#!/bin/sh

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d    : default"
    echo " [0] -h    : print help"
}

if [ $# -eq 0 ]
then
    printHelp
else
    if [ "$1" = "-d" ]; then
        rnd_seed=`date +%N`
	nev=100000000
	statisticsMultiplyFactor=10
	cmd="./cosmique_gamma_hadron_generator 0 $nev cosmique_gamma_hadron_generator_h200km.root $rnd_seed $statisticsMultiplyFactor"
        #echo $cmd
        $cmd
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi
