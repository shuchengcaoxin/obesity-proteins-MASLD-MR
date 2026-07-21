#!/bin/bash
#PBS -N MR #Name for the job
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID} #path for output
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID} #path for error
#PBS -l walltime=12:00:00 #Maximum time for the job
#PBS -l nodes=1:ppn=1 #specified number of nodes with a certain number of processors per node
#PBS -l mem=24G #Maximum amount of RAM
#PBS -l vmem=24G #Virtual memory
#PBS -t 1-491 #specifies a job array with tasks numbered from 1 to 491
 
mkdir -p output/harmonized #The '-p' option tells 'mkdir' to make parent diretories as necessary
mkdir -p output/or
mkdir -p output/pleio
mkdir -p output/hetero

cd $PBS_O_WORKDIR #The absolute path of the current working directory of the qsub command
mkdir -p logs/

list=$(head -n ${PBS_ARRAYID} /scratch/richards/satoshi.yoshiji/11.pQTL/decode_batch/listbatch.txt | tail -n 1) # path to pQTLs list
#head -n ${PBS_ARRAYID}, ID is a number (N), it prints N lines of the file
#/scratch/.../.txt the file to be read
#tail -n 1 return the last line of that output from 'head' command
#so every time it gives a specific Nth line of the file and put it into the list

for protein in `cat $list`; #cat: display the contents of a file on the screen
do
protein_name=`basename $protein ".txt.gz"`#It is a bit hard for me to understand if I didn't see the file

/admin/kleinman/custom_software/R-4.1.2-gcc-9.3.0/lib64/R/bin/Rscript --vanilla 01.MR.R $protein $protein_name
#--vanilla tells 'Rscript' to start R with a minimal set of settings
#Avoid loading site-specific initialization scripts, no saed workspae, no command history.
#Ensure the sript runs in a clean and reproducible environment, without previous set options' interference
done
