from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

gro = GromacsGroFile('solv_ions.gro')
top = GromacsTopFile('complex.top', periodicBoxVectors=gro.getPeriodicBoxVectors(), includeDir='/usr/local/gromacs/share/gromacs/top')

system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

platform = Platform.getPlatformByName('CUDA')
properties = {'CudaDeviceIndex':'0', 'CudaPrecision':'mixed'}


simulation = Simulation(top.topology, system, integrator, platform, properties)
simulation.context.setPositions(gro.positions)
simulation.minimizeEnergy()
simulation.reporters.append(DCDReporter('output-2ndgpu-mp.dcd', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, progress=True, remainingTime=True, speed=True, totalSteps=50000, separator='\t'))
simulation.step(50000)
