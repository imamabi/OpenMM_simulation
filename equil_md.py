"""
Run a MD simulation for a complex, optionally adding a solvent box

This script was obtained from https://github.com/tdudgeon/simple-simulate-complex

Few adjustments were made to perform NPT equilibration and NVT production

"""

import sys, time, argparse

from openff.toolkit import Molecule
from openmmforcefields.generators import SystemGenerator
import openmm
from openmm import *
from openmm import app, unit, LangevinIntegrator, Vec3
from openmm.app import PDBFile, Simulation, Modeller, PDBReporter, StateDataReporter, DCDReporter, CheckpointReporter

import utils

t0 = time.time()


parser = argparse.ArgumentParser(description="simulateComplex")

parser.add_argument("-p", "--protein", required=True, help="Protein PDB file")
parser.add_argument("-l", "--ligand", required=True, help="Ligand molfile")
parser.add_argument("-o", "--output", default='output', help="Base name for output files")
parser.add_argument("-z", "--step-size", type=float, default=0.002, help="Step size (ps")
parser.add_argument("-f", "--friction-coeff", type=float, default=1, help="Friction coefficient (ps)")
parser.add_argument("-i", "--interval", type=int, default=100000, help="Reporting interval")
parser.add_argument("-t", "--temperature", type=int, default=300, help="Temperature (K)")
parser.add_argument("--solvate", action='store_true', help="Add solvent box")
parser.add_argument("--padding", type=float, default=15, help="Padding for solvent box (A)")
parser.add_argument("--water-model", default="tip3p",
                    choices=["tip3p", "spce", "tip4pew", "tip5p", "swm4ndp"],
                    help="Water model for solvation")
parser.add_argument("--positive-ion", default="Na+", help="Positive ion for solvation")
parser.add_argument("--negative-ion", default="Cl-", help="Negative ion for solvation")
parser.add_argument("--ionic-strength", type=float, default="0.15", help="Ionic strength for solvation")
parser.add_argument("--no-neutralize", action='store_true', help="Don't add ions to neutralize")
parser.add_argument("-s", "--equilibration-steps", type=int, default=200, help="Number of equilibration steps")
parser.add_argument("--protein-force-field", default='amber/ff14SB.xml', help="Metallic ion force field")
parser.add_argument("--metal-ion-force-field", default='amber/tip3p_HFE_multivalent.xml', help="Protein force field")
parser.add_argument("--ligand-force-field", default='gaff-2.11', help="Ligand force field")
parser.add_argument("--water-force-field", default='amber/tip3p_standard.xml', help="Ligand force field")

args = parser.parse_args()
print("simulateComplex: ", args)

pdb_in = args.protein
mol_in = args.ligand
output_base = args.output
output_complex = output_base + '_complex.pdb'
output_traj_dcd = output_base + '_traj.dcd'
output_min = output_base + '_minimised.pdb'
reporting_interval = args.interval
temperature = args.temperature * unit.kelvin
equilibration_steps = args.equilibration_steps
print('Processing', pdb_in, 'and', mol_in, 'with', equilibration_steps, 'steps generating outputs',
      output_complex, output_min, output_traj_dcd)

# get the chosen or fastest platform
platform = utils.get_platform()

print('Reading ligand')
ligand_mol = Molecule.from_file(mol_in)

print('Preparing system')
# Initialize a SystemGenerator using the GAFF for the ligand and tip3p for the water.
forcefield_kwargs = {'constraints': app.HBonds, 'rigidWater': True, 'removeCMMotion': False, 'hydrogenMass': 4*unit.amu }
system_generator = SystemGenerator(
    forcefields=[args.protein_force_field, args.water_force_field, args.metal_ion_force_field],
    small_molecule_forcefield=args.ligand_force_field,
    molecules=[ligand_mol],
    forcefield_kwargs=forcefield_kwargs)

# Use Modeller to combine the protein and ligand into a complex
print('Reading protein')
protein_pdb = PDBFile(pdb_in)

print('Preparing complex')
# The topology is described in the openforcefield API

modeller = Modeller(protein_pdb.topology, protein_pdb.positions)
print('System has %d atoms' % modeller.topology.getNumAtoms())

# The topology is described in the openforcefield API
print('Adding ligand...')
lig_top = ligand_mol.to_topology()
modeller.add(lig_top.to_openmm(), lig_top.get_positions().to_openmm())
print('System has %d atoms' % modeller.topology.getNumAtoms())

# Solvate
if args.solvate:
    print('Adding solvent...')
    # we use the 'padding' option to define the periodic box.
    # we just create a box that has a 10A (default) padding around the complex.
    modeller.addSolvent(system_generator.forcefield, model=args.water_model, padding=args.padding * unit.angstroms,
                        positiveIon=args.positive_ion, negativeIon=args.negative_ion,
                        ionicStrength=args.ionic_strength * unit.molar, neutralize=not args.no_neutralize)
    print('System has %d atoms' % modeller.topology.getNumAtoms())

with open(output_complex, 'w') as outfile:
    PDBFile.writeFile(modeller.topology, modeller.positions, outfile)

# Create the system using the SystemGenerator
system = system_generator.create_system(modeller.topology, molecules=ligand_mol)

friction_coeff = args.friction_coeff / unit.picosecond
step_size = args.step_size * unit.picoseconds
duration = (step_size * equilibration_steps).value_in_unit(unit.nanoseconds)
print('Equilibrating for {} ns'.format(duration))

integrator = LangevinIntegrator(temperature, friction_coeff, step_size)
if args.solvate:
    system.addForce(openmm.MonteCarloBarostat(1 * unit.atmospheres, temperature, 25))

if system.usesPeriodicBoundaryConditions():
    print('Default Periodic box: {}'.format(system.getDefaultPeriodicBoxVectors()))
else:
    print('No Periodic Box')

simulation = Simulation(modeller.topology, system, integrator, platform=platform)
context = simulation.context
context.setPositions(modeller.positions)

print('Minimising ...')
simulation.minimizeEnergy()

# Write out the minimised PDB.
with open(output_min, 'w') as outfile:
    PDBFile.writeFile(modeller.topology, context.getState(getPositions=True, enforcePeriodicBox=True).getPositions(), file=outfile, keepIds=True)

# Write out the system information.
with open('system.xml', 'w') as output:
    output.write(XmlSerializer.serialize(system))
    
# equilibrate
simulation.context.setVelocitiesToTemperature(temperature)
simulation.reporters.append(StateDataReporter(sys.stdout, reporting_interval * 5, step=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(CheckpointReporter('equilibration_checkpnt.chk', (reporting_interval/20)))
simulation.reporters.append(StateDataReporter("equilibration_log.csv", reporting_interval, step=True, potentialEnergy=True, temperature=True, volume=True))

print('Equilibrating ...')
simulation.step(equilibration_steps)

# Save topology as a pdb file
with open('topology.pdb', 'w') as output:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(),file=output, keepIds=True)

t1 = time.time()
print('NPT equilibration complete in {} mins at {}.'.format(
    round((t1 - t0) / 60, 3), temperature))
print('Equilibration time was', round(duration, 3), 'ns')
