import openmm
from openmm import *
from openmm import app, unit, LangevinIntegrator, Vec3
from openmm.app import PDBFile, Simulation, Modeller, PDBReporter, StateDataReporter, DCDReporter, CheckpointReporter

import sys, time, argparse

import utils

parser = argparse.ArgumentParser(description="simulateComplex")

parser.add_argument("-o", "--output", default='output', help="Base name for output files")
parser.add_argument("-s", "--steps", type=int, default=5000, help="Number of steps")
parser.add_argument("-z", "--step-size", type=float, default=0.002, help="Step size (ps")
parser.add_argument("-f", "--friction-coeff", type=float, default=1, help="Friction coefficient (ps)")
parser.add_argument("-i", "--interval", type=int, default=100000, help="Reporting interval")
parser.add_argument("-t", "--temperature", type=int, default=300, help="Temperature (K)")

args = parser.parse_args()
print("simulateComplex: ", args)

output_base = args.output
output_traj_dcd = output_base + '_traj.dcd'
num_steps = args.steps
reporting_interval = args.interval
temperature = args.temperature * unit.kelvin

friction_coeff = args.friction_coeff / unit.picosecond
step_size = args.step_size * unit.picoseconds

platform = utils.get_platform()

# Building your simulation system
#import saved files
model_sys = PDBFile('topology.pdb')

with open('system.xml') as input:
    system = XmlSerializer.deserialize(input.read())

# Define integrators
#integrator = LangevinMiddleIntegrator(temperature, friction_coeff, step_size)
integrator = VerletIntegrator(step_size)
system.addForce(openmm.AndersenThermostat(temperature, friction_coeff))

#create simulation
simulation=Simulation(model_sys.topology, system, integrator, platform=platform)

# Function to load a checkpoint and return the simulation context
def load_checkpoint(checkpoint_file, simulation):
    with open(checkpoint_file, 'rb') as file:
        simulation.context.loadCheckpoint(file.read())
    return simulation.context

# Function to calculate remaining steps needed
def calculate_remaining_steps(total_steps, simulation):
    steps_completed = simulation.context.getState().getStepCount()
    return max(0, total_steps - steps_completed)

# Load last saved checkpoint file
checkpoint_file = 'production_checkpnt.chk'
simulation_context = load_checkpoint(checkpoint_file, simulation)
# Total desired number of steps for the simulation
total_steps_desired = num_steps

# Calculate remaining steps needed
remaining_steps = calculate_remaining_steps(total_steps_desired, simulation)

# Get duration for simulation
duration = (step_size * remaining_steps).value_in_unit(unit.nanoseconds)

#Create reporters
simulation.reporters.append(DCDReporter(output_traj_dcd, reporting_interval, enforcePeriodicBox=True, append=True))
simulation.reporters.append(StateDataReporter(sys.stdout, reporting_interval * 5, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(CheckpointReporter('production_checkpnt.chk', (reporting_interval/20)))
simulation.reporters.append(StateDataReporter("production_log.csv", reporting_interval, step=True, potentialEnergy=True, temperature=True, volume=True, append=True))

print('Continue NVT simulation to', num_steps, 'steps ...')
t1 = time.time()
simulation.step(remaining_steps)
t2 = time.time()

print('Simulation complete in {} mins at {}.'.format(
    round((t2 - t1) / 60, 3), temperature))
print('Simulation time was', round(duration, 3), 'ns')
