#!/bin/sh

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d    : setup environment"
    echo " [0] -l    : load modules"
    echo " [0] -p    : purge modules"
    echo " [0] -list : list modules"
    echo " [0] -h    : print help"
}

if [ $# -eq 0 ] 
then    
    printHelp
else
    if [ "$1" = "-d" ]; then
	module load GCC/8.3.0 OpenMPI/3.1.4 ROOT/6.20.04-Python-3.7.4
	source /opt/ebsofts/ROOT/6.20.04-foss-2019b-Python-3.7.4/bin/thisroot.sh
    elif [ "$1" = "-l" ]; then
	module load GCC/8.3.0 OpenMPI/3.1.4 ROOT/6.20.04-Python-3.7.4
    elif [ "$1" = "-p" ]; then
	module purge
    elif [ "$1" = "-list" ]; then
	echo " "
	echo " "
	echo " "
	echo " "
	echo "module list"
	module list
	echo " "
	echo " "
	echo " "
	echo " "
	echo "root"
	printenv | grep -i root | head -n 5 | grep -i root
	#printenv | grep -i root
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi

#espeak "I have done"
