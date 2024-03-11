from openmm.app import *
from openmm import *
from openmm.unit import *
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda
from sys import stdout
import os

path_to_sim = os.getcwd()
os.chdir(path=path_to_sim)  # Set the working directory to the simulation files directory

#Choose simulation platform and set precision
platform = utils.get_platform()

# Import amber coordinate and topology files
print('Loading...')

files = os.listdir(path_to_sim)
for file in files:
    if file.endswith(".inpcrd"):
        inpcrd = AmberInpcrdFile(f'{path_to_sim}/{file}')
    
    elif file.endswith(".prmtop"):
        prmtop = AmberPrmtopFile(f'{path_to_sim}/{file}', periodicBoxVectors=inpcrd.boxVectors)
    
    else:
        continue

# Create system
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)

# Write out the system information.
with open('system.xml', 'w') as output:
    output.write(XmlSerializer.serialize(system))

# Minimize energy
print('Minimizing...')
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds, platform=platform)

#Create simulation
simulation = Simulation(prmtop.topology, system, integrator, platform=platform)
simulation.context.setPositions(inpcrd.positions)

#Carry out energy minization
try:
    simulation.minimizeEnergy(maxIterations=50000, tolerance=10)
    print('Minimization completed successfully.')
except Exception as e:
    print(f'Minimization failed: {e}')

# Write out the minimised PDB.
with open('system_minimized.pdb', 'w') as outfile:
    PDBFile.writeFile(simulation.topology, context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(), 
                      file=outfile, keepIds=True)

# NVT equilibration
print('NVT equilibration starting...')
equilibration_reporters = [StateDataReporter("equil_nvt_log.txt", 5000, step=True, 
                                             potentialEnergy=True, temperature=True, volume=True),
                                             CheckpointReporter('equil_nvt.chk', 10000)]

simulation.reporters = equilibration_reporters
simulation.step(10000000)

# NPT equilibration
print('NPT equilibration starting...')
system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
simulation.context.reinitialize(preserveState=True)

equil_reporters = [StateDataReporter("equil_npt_log.txt", 5000, step=True, 
                                             potentialEnergy=True, temperature=True, volume=True),
                                             CheckpointReporter('equil_npt.chk', 10000)]

simulation.reporters = equil_reporters
simulation.step(10000000)

#import saved files
model_sys = PDBFile('system_minimized.pdb')

with open('system.xml') as input:
    system = XmlSerializer.deserialize(input.read())

#define integrators
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

#create simulation
simulation=Simulation(model_sys.topology, system, integrator, platform=platform)

#set positions, velocities and box vectors from checkpoint
simulation.loadCheckpoint('equil_npt.chk')

# Production MD
print('Production MD starting...')
production_reporters = [DCDReporter('production_traj.dcd', 100000), CheckpointReporter('production_checkpnt.chk', 10000)
                        StateDataReporter("production_log.txt", 5000, step=True,
                                          potentialEnergy=True, temperature=True, volume=True)]
simulation.reporters = production_reporters
simulation.step(50000000)  # Adjust the number of steps as needed

print('Simulation finished!')