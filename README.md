# Argon-Simulator

These instructions provide details on how to run the code.

Modules

argon_simulator.py – contains the main method to run the simulation as well as functions to produce MSD and RDF data.

Particle3D.py – a class to create the particles and update properties.

mdutilities_20210131.py – a module to initialise initial parameters (lattice size, positions, velocities).

pbc.py – module containing functions to return particle positions in periodic boundary conditions. 

ljpotenyial.py – module to apply Lennard Jones potential condition.

Running Test (in terminal)

To run the test for e.g. a solid:

python argon_simulator.py solid_traj.xyz solid_obs.txt particles_solid.txt param_solid.txt

To run for a liquid and gas:

python argon_simulator.py liquid_traj.xyz liquid_obs.txt particles_liquid.txt param_liquid.txt

python argon_simulator.py gas_traj.xyz gas_obs.txt particles_gas.txt param_gas.txt

These will produce the following files:

[state]_traj.xyz – contains position data at each timestep to load in to VMD for visualisation.

[state]_obs – contains total energy, kinetic energy, potential energy, and MSD data for each timestep.

State being the solid, liquid or gas.

The program will also create plots for MSD, RDF, and the energies. 

Units used are reduced units: distance (σ), energy (ϵ), mass (m_0), time (τ=σ √((m/ϵ) )).
