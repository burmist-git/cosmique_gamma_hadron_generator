#!/bin/sh

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d         : default"
    echo " [0] -c         : corsika"
    echo " [0] -con       : corsika on axis"
    echo " [0] -recompile : recompile"
    echo " [0] -h         : print help"
}

if [ $# -eq 0 ]
then
    printHelp
else
    if [ "$1" = "-d" ]; then
        rnd_seed=`date +%N`
	nev=10000
	statisticsMultiplyFactor=1
	cmd="./cosmique_gamma_hadron_generator 0 $nev cosmique_gamma_hadron_generator.root $rnd_seed $statisticsMultiplyFactor"
        #echo $cmd
        $cmd
    elif [ "$1" = "-c" ]; then
        rnd_seed=`date +%N`
	nev=100000000
	statisticsMultiplyFactor=1
	cmd="./cosmique_gamma_hadron_generator 0 $nev cosmique_gamma_hadron_generator_corsika.root $rnd_seed $statisticsMultiplyFactor"
        $cmd
    elif [ "$1" = "-con" ]; then
        rnd_seed=`date +%N`
	nev=100000000
	statisticsMultiplyFactor=1
	cmd="./cosmique_gamma_hadron_generator 0 $nev cosmique_gamma_hadron_generator_corsika_onaxis.root $rnd_seed $statisticsMultiplyFactor"
        $cmd
    elif [ "$1" = "-recompile" ]; then
	make clean ; make;
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi
