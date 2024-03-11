#!/bin/bash


#Remember to install openmm, openmmforcefields and its dependencies.


/Path/to/conda activate openmm

pyhton prepareProtein protein.pdb output

#Execute the equilibration script
python equil_md.py -p output_fixed.pdb -l ligand.sdf -o prot_lig -s 100000  -z 0.002 -i 100000 -t 310 --solvate

#Execute the production script
python prod_md.py -o prot_lig.pdb -s 250000  -z 0.002 -i 100000 -t 310


