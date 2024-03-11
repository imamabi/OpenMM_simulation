import sys, time, argparse

from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator
import openmm
from openmm import *
from openmm import app, unit, LangevinIntegrator, Vec3
from openmm.app import PDBFile, Simulation, Modeller, PDBReporter, StateDataReporter, DCDReporter, CheckpointReporter

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

duration = (step_size * num_steps).value_in_unit(unit.nanoseconds)

# Building your simulation system
#import saved files
model_sys = PDBFile('topology.pdb')

with open('system.xml') as input:
    system = XmlSerializer.deserialize(input.read())

#define integrators
integrator = LangevinMiddleIntegrator(temperature, friction_coeff, step_size)

#create simulation
simulation=Simulation(model_sys.topology, system, integrator, platform=platform)

#set positions, velocities and box vectors from checkpoint
simulation.loadCheckpoint('equilibration_checkpnt.chk')

#Create reporters
simulation.reporters.append(DCDReporter(output_traj_dcd, reporting_interval, enforcePeriodicBox=True))
simulation.reporters.append(StateDataReporter(sys.stdout, reporting_interval * 5, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(CheckpointReporter('production_checkpnt.chk', (reporting_interval/20)))
simulation.reporters.append(StateDataReporter("prodcution_log.csv", reporting_interval, step=True, potentialEnergy=True, temperature=True, volume=True))

print('Starting NVT simulation with', num_steps, 'steps ...')
t1 = time.time()
simulation.step(num_steps)
t2 = time.time()

print('Simulation complete in {} mins at {}.'.format(
    round((t2 - t1) / 60, 3), temperature))
print('Simulation time was', round(duration, 3), 'ns')
