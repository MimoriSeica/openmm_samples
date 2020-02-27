from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

def annealing(integrator, simulation, target_temperature):
    temperature_list = list(range(0, target_temperature, 3))
    temperature_list.append(target_temperature)
    for temperature in temperature_list:
        integrator.setTemperature(temperature*kelvin)
        simulation.step(100)


pdb = PDBFile('model/alat.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()

"""
annealing目標の温度まで三度つづ行う
minimizeEnergy() の後に行う
"""
annealing(integrator, simulation, 300)

simulation.reporters.append(PDBReporter('output.pdb', 100))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000)
