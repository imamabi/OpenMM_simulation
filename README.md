# OpenMM_simulation
## Overview
This project focuses on the simulation of protein-ligand complexes in explicit solvent using OpenMM. The workflow includes NPT equilibration, NVT production, and extending simulation time, providing comprehensive insights into the molecular dynamics of the system.

## Workflow
### 1. NPT Equilibration
Script: equil_md.py

Description: Prepares the protein-ligand complex (assumes the ligand is properly docked to the receptor) and performs NPT (constant Number of particles, Pressure, and Temperature) equilibration.


### 2. NVT Production
Script: prod_md.py

Description: Continues from the last equilibration checkpoint and performs NVT (constant Number of particles, Volume, and Temperature) molecular dynamics production.


### 3. Extension of Simulation Time
Script: cont_md.py

Description: Extends the simulation time from the previous production run.

## Repository Contents
MDrun-openmm-lcc.sh: A bash script that automates the simulation process from energy minimization (EM) to production on a high-performance computing (HPC) cluster.

README.md: Project documentation.

cont_md.py: Script to extend the simulation time.

equil_md.py: Script for NPT equilibration of the protein-ligand complex.

openMM_amber.py: Script for setting up the simulation using Amber force fields.

openMM_bashrun.sh: A bash script to run OpenMM simulations.

prepareProtein.py: A PDBFixer script that preprocesses the protein PDB file.

prod_md.py: Script for NVT molecular dynamics production.

utils.py: Utility functions used across different scripts.
