#!/bin/bash

# Create Working Directory
WDIR=$SCRATCH/$JOB_NAME-$JOB_ID-$3
ROSETA=$WORK/MyProjects/rosetta


mkdir -p $WDIR
if [ ! -d $WDIR ]
then
  echo $WDIR not created
  exit
fi
cd $WDIR

# Copy Data and Config Files
cp $WORK/MyProjects/UCSF_RSA_runs/config/* .
cp $WORK/MyProjects/UCSF_RSA_runs/structs_art/$1.pdb .
# Put your Science related commands here
date
hostname


date
#This is the line that runs fixbb method that makes the mutations
./fixbb -database /home1/02212/eleishaj/ROSETTA_DB/ -s $1.pdb -resfile ALLAA.res -ex1 -ex2 -extrachi_cutoff 0 -nstruct 1 -out::prefix $1_$3 -overwrite -linmem_ig 10 >& fixbb.log.txt
date

# Copy Results Back to Home Directory
RDIR=$HOME/MyProjects/UCSF_RSA_runs/$1/yeast_RSA-$1-$JOB_ID-$3
mkdir -p $RDIR
cp * $RDIR/. 

# Cleanup 
rm -rf $WDIR
