#!/bin/bash

#SBATCH -t 72:00:00                                   #Time for the job to run
#SBATCH --job-name=G12C_GDP.                          #Name of the job
#SBATCH -N 1                                          #Number of nodes required
#SBATCH --ntasks 10                                   #Number of cores needed for the job
##SBATCH --partition=V4V16_SKY32M192_L                 #Name of the queue
#SBATCH --partition=V4V32_CAS40M192_L                           #Name of the queue
##SBATCH --partition=V4V32_SKY32M192_L                           #Name of the queue
#SBATCH --gres=gpu:1                                  #Number of GPU's
#SBATCH --mail-type ALL                               #Send email on start/end
##SBATCH --nodelist=gvnodeb004
#SBATCH --mail-type ALL                               #Send email on start/end
#SBATCH --mail-user iaim223@uky.edu                   #Where to send email
#SBATCH --account=gol_qsh226_uksr                     #Name of account to run under


module --force purge
module load gnu12/12.1.0
module load ccs/cuda/10.0.130 
#module load ccs/cuda/11.2.0_460.27.04
#module load ccs/cuda/11.6.0_510.39.01

source activate /project/qsh226_uksr/iaim223/openmm_env


#Execute the equilibration script
python equil_md.py -p G12C_fixed.pdb -l gdp.sdf -o G12C_GDP -s 100000000  -z 0.002 -i 100000 -t 310 --solvate

#Execute the production script
python prod_md.py -o G12C_GDP -s 500000000  -z 0.002 -i 100000 -t 310

# Extend simulation from last checkpoint
#python cont_md.py -o G12C_GDP -s 1000000000  -z 0.002 -i 100000 -t 310              
