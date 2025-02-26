# OpenMM_simulation
## Overview
This project focuses on the simulation of protein-ligand complexes in explicit solvent using OpenMM. The workflow includes NPT equilibration, NVT production, and extending simulation time, providing comprehensive insights into the molecular dynamics of the system.

## Workflow
### 1. NPT Equilibration
Script: equil_md.py

Description: Prepares the protein-ligand complex (assumes the ligand is properly docked to the receptor) and performs NPT (constant Number of particles, Pressure, and Temperature) equilibration.

Output: Equilibrated structures of the protein-ligand complex.

### 2. NVT Production
Script: prod_md.py

Description: Continues from the last equilibration checkpoint and performs NVT (constant Number of particles, Volume, and Temperature) molecular dynamics production.

Output: Production trajectories of the protein-ligand complex.
