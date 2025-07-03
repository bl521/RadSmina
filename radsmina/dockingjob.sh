#!/bin/bash

#PBS -lwalltime=28:00:00 
#PBS -lselect=1:ncpus=8:mem=64gb 

source ~/anaconda3/bin/activate radsmina
 
cd $PBS_O_WORKDIR

python3 $PBS_O_WORKDIR/DUDEZ_smina.py