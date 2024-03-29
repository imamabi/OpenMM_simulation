# OpenMM_simulation
Simulation of protein ligand complex in explicit solvent (NPT equilibration and NVT production)
The equil_md.py prepares the protein ligand complex (assumes ligand is properly docked to the receptor). and performs NPT equilibration.
The prod_md.py continues from the last equilibration checkpoint and performs NVT MD production.
The cont_md.py is use to extent the simulation time.
