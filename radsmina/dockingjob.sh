#!/bin/bash

#PBS -lwalltime=28:00:00 
#PBS -lselect=2:ncpus=64:mem=64gb 

source ~/anaconda3/bin/activate rad1
 
cd /rds/general/user/bl521/home/rad/examples

python3 smina_dock.py