#!/bin/sh
#SBATCH --job-name easc%j
#SBATCH --error /srv/beegfs/scratch/users/b/burmistr/cosmique_gamma_hadron_generator/job_error/crgen_%j.error
#SBATCH --output /srv/beegfs/scratch/users/b/burmistr/cosmique_gamma_hadron_generator/job_output/output_%j.output
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 1
#SBATCH --partition public-cpu
#SBATCH --time 10:00:00

function printHelp {
    echo " --> ERROR in input arguments "
    echo " [0] -d    : single job"
    echo " [1]       : jobID"
    echo " [0] -h    : print help"
}

if [ $# -eq 0 ] 
then
    printHelp
else
    if [ "$1" = "-d" ]; then
        if [ $# -eq 2 ]; then
            jobID=$2
            job_out="/srv/beegfs/scratch/users/b/burmistr/cosmique_gamma_hadron_generator/out/$jobID"
            mkdir -p $job_out
            rnd_seed=`date +%N`
	    nev=1000000000
	    statisticsMultiplyFactor=10
	    cmd="/home/users/b/burmistr/cosmique_gamma_hadron_generator/cosmique_gamma_hadron_generator 0 $nev $job_out/cosmique_gamma_hadron_generator.root $rnd_seed $statisticsMultiplyFactor"
            #echo $cmd
            $cmd
        else
            printHelp
        fi
    elif [ "$1" = "-h" ]; then
        printHelp
    else
        printHelp
    fi
fi


