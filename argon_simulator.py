# -*- coding: utf-8 -*-
"""
Argon Simulator; code to create a system of Argon particles and caluclate its 
energy, MSD and RDF. 

Units used are reduced units: distance (σ), energy(ϵ), mass (m_0), time 
(τ=σ √((m/ϵ) )).

@author: Isaac Watson
"""

import numpy as np
from Particle3D import Particle3D
import sys
import ljpotential
import mdutilities_20210131 as mdutilities
import pbc
import matplotlib.pyplot as pyplot

def bin_data(positions, bins, l):
    """
    Retrieves data to put in the histogram bins
    
    :param positions : current positions (array)
    :param bins : number of bins (integer)
    :param l : lattice size (vector)
    :return particle_count : raw particle count data (integer)
    :return bin_edges : edges of the bins (array)

    """
    N = len(positions)
    r_ij = []
    #Compute separations
    for i in range(N):
        for j in range(i):
            if i!= j:
                r_mic = Particle3D.mic_sep(positions[i], positions[j], l)
                r_ij.append(np.linalg.norm(r_mic))
    # Get arrays of the values of the historgam and the bin boundaries
    # Divide by N/2 as part of normalisation and because only half separations 
    # considered(i-j not j-i)
    particle_count = np.histogram(r_ij, bins)[0]/(N/2)
    bin_edges = np.histogram(r_ij, bins)[1]
    
    return particle_count, bin_edges
    
def rdf_norm(positions, rho, bins, l, count_avg):
    """
    Creates data to produce a normalised histogram 
    
    :param positions : current positions (array)
    :param rho : desnity (float)
    :param bins : number of bins (integer)
    :param l : lattice size (vector)
    :param count_avg : time average of count 
    :return rdf_hist : histogram values (float)
    :return bin_mids : histogram bin mid points (array)
    """
    # Get histogram data and bin edges
    hist_data, bin_edges = bin_data(positions, bins, l)
    # Find midpoints of the bins (r)
    bin_mids = 0.5*(bin_edges[1:]+bin_edges[:-1])
    # Find bin width (delta_r)
    bin_width = bin_edges[1]-bin_edges[0]
    # Calculate the expected historgram values
    expected = 4*np.pi*rho*bin_mids**2*bin_width
    # Calculate the rdf (g(r)) 
    rdf_hist = count_avg/(expected) 
    
    return rdf_hist, bin_mids
    
    
def msd(init_pos, positions, l):
    """
    Retrieves MSD data for the particles 
    
    :param init_pos : inital positions (array)
    :param positions : current positions (array)
    :param l : lattice size (vector)
    :return msd : mean squared displacement (float)

    """
    # Create empty displacment array
    delta_r = np.zeros([len(init_pos), 3])
    # For each element in the array calculate the difference 
    # between the initial and current positions (displacement)
    for i in range(len(init_pos)):
        delta_r[i] = pbc.mic(positions[i]-init_pos[i], l)
    # Sum the squares of the displacemnts
    sum_sq = np.sum(delta_r**2)
    # Divide by number of particles to calculate MSD
    msd = sum_sq/(len(init_pos))
    
    return msd
    

def vmd(xyz, p3d_list, i):
    """
    Writes xyz position data to an xyz trajectory file for VMD

    :param xyz: file name (string)
    :param p3d_list: list of particles (list)
    :param i: numstep value (integer)
        
    """
    
    # Write frist two lines to match xyz file 
    xyz.write(str(len(p3d_list)) + "\n")
    xyz.write(str(i) + "\n")
    # Write xyz string for all particles in particle list
    for p in p3d_list:
        xyz.write(str(p) + "\n")


# Begin main code
def main():
    # Read name of output files from command line
    #if len(sys.argv)!=5:
        #print("Wrong number of arguments.")
        #print("Usage: " + sys.argv[0] + " <output file>")
        #quit()
    #else:
    outfile_name = "melting_traj1.21.xyz" #sys.argv[1]
    outfile2_name = "melting_obs1.21.txt" #sys.argv[2]
    particle_data = "particles_solid.txt" #sys.argv[3]
    param_data = "param_melting.txt" #sys.argv[4]

        

    # Open output file
    outfile = open(outfile_name, "w")
    outfile2 = open(outfile2_name, "w")


    # Open parameter data file
    params = open(param_data, "r")
    # Make sure file is in read mode
    if params.mode == "r":
        # Read first line and split into a list so items can be assigned a variable
        line = params.readline().split(" ")
        # Assign variables to items to set up simulation parameters 
        dt = float(line[0])
        numstep = int(line[1])
        temp = float(line[2])
        rho = float(line[3])
        cutoff = float(line[4])
    
    # Set up particle initial conditions 
    particles = open(particle_data, "r")
    # Make sure file is in read mode
    if particles.mode == "r":
        # Read first line
        line = particles.readline().split(" ")
        # Read in number of partcicles
        N = int(line[0])
    
    # Create list of partcicles
    p3d_list = []
    for i in range(1, N+1):
        p3d_list.append(Particle3D.new_particle(i))
     
    # Rertieve lattice size, initial positions and initial velocities    
    l = mdutilities.set_initial_positions(rho, p3d_list)[0] 
    mdutilities.set_initial_velocities(temp, p3d_list)

    # Create an array of intial positions to keep fixed
    init_positions = np.zeros([N, 3])
    for i in range(N):
        init_pos = p3d_list[i].p
        init_positions[i] = init_pos

    # Create array of initial positions to update
    positions = np.zeros([N, 3])
    for i in range(N):
        pos = p3d_list[i].p
        positions[i] = pos

    # Create array of intial forces
    forces_array = ljpotential.get_forces(init_positions, l, cutoff)
    
    # Determine bins for rdf
    bins = np.arange(0, int(l[0]), 0.1)
    
    # Write out initial conditions
    msd_val = msd(init_positions, positions, l)
    k_energy = Particle3D.sys_kinetic(p3d_list)
    p_energy = ljpotential.sys_potential_lj(positions, l, cutoff) 
    t_energy = k_energy + p_energy
    time = 0.0
    
    
    # Create lists to append data to for graphs
    energy_list = [t_energy]
    pe_list = [p_energy]
    ke_list = [k_energy]
    time_list = [time]
    msd_list = [msd_val]

    count, bins = bin_data(positions, bins, l)
 
    # Write observables to file
    outfile2.write("{} {} {} {} {} \n".format(time, t_energy, k_energy, p_energy, msd_val))
    
    # Start the time integration loop
    for i in range(numstep):
        
        # Print i at regular intervals to monitor codes progress
        if i % 100 == 0:
            print(i, numstep)
        
        # Update particle positions (second order update)
        Particle3D.update_pos_all(p3d_list, dt, forces_array, l)
        
        # Assign new positions to array
        for j in range(N):
            pos = p3d_list[j].p          
            positions[j] = pos

        # Update forces (create new array)
        forces_new = ljpotential.get_forces(positions, l, cutoff)
        
        # Update particle velocity by averaging current and new forces
        Particle3D.update_vel_all(p3d_list, dt, 0.5*(forces_array+forces_new))

        # Re-define force values in force array
        forces_array = forces_new
        
        # Increase time
        time += dt
            
        # Calculate energies
        k_energy = Particle3D.sys_kinetic(p3d_list)
        p_energy = ljpotential.sys_potential_lj(positions, l, cutoff) 
        t_energy = k_energy + p_energy
        
        # Get data for MSD
        msd_val = msd(init_positions, positions, l)
        
        #if i % 100 == 0:
        #Retrieve particle count data for RDF
        new_count = bin_data(positions, bins, l)[0]
        count += new_count

        # Append data values to lists
        energy_list.append(t_energy)
        pe_list.append(p_energy)
        ke_list.append(k_energy)
        time_list.append(time)
        msd_list.append(msd_val)

        # Write updated observables to file
        outfile2.write("{} {} {} {} {} \n".format(time, t_energy, k_energy, p_energy, msd_val))
        
        # Write VMD xyz file
        vmd(outfile, p3d_list, i)
     
    # Calculate time average of particle count data
    count_avg = count/numstep 
    # Retrieve data fro graphing
    rdf_data, bin_points = rdf_norm(positions, rho, bins, l, count_avg)
    # Post-simulation:
    # Close output files
    outfile.close()
    outfile2.close()
    

    # Plot all particle energies to screen
    pyplot.title('Argon: energy vs time')
    pyplot.xlabel('Time/tau')
    pyplot.ylabel('Energy/epsilon')
    pyplot.plot(time_list, energy_list, color = 'red')
    pyplot.plot(time_list, pe_list, color = 'blue')
    pyplot.plot(time_list, ke_list, color = 'green')
    pyplot.legend(['Total', 'Potential', 'Kinetic'])
    pyplot.show()
    
    # Plot total particle energy to screen
    pyplot.title('Argon: total energy vs time')
    pyplot.xlabel('Time/tau')
    pyplot.ylabel('Energy/epsilon')
    pyplot.plot(time_list, energy_list, color = 'red')
    pyplot.show()
    
    # Plot KE
    pyplot.title('Argon: kinetic energy vs time')
    pyplot.xlabel('Time/tau')
    pyplot.ylabel('Energy/epsilon')
    pyplot.plot(time_list, ke_list, color = 'green')
    pyplot.show()
    
    # Plot PE
    pyplot.title('Argon: potential energy vs time')
    pyplot.xlabel('Time/tau')
    pyplot.ylabel('Energy/epsilon')
    pyplot.plot(time_list, pe_list, color = 'blue')
    pyplot.show()
    
    # Plot MSD to screen   
    pyplot.title('Argon: msd vs time')
    pyplot.xlabel('Time/tau')
    pyplot.ylabel('MSD/sigma^2')
    pyplot.plot(time_list, msd_list)
    pyplot.show()
    
    # Plot RDF to screen
    pyplot.title('Argon: rdf vs r')
    pyplot.xlabel('r/sigma')
    pyplot.ylabel('rdf')
    # Sum RDF data in each bin over all numsteps
    pyplot.plot(bin_points, rdf_data)
    pyplot.show()
    

    outfile3 = open("melting_rdf1.21.txt", "w")
    c = [bin_points, rdf_data]
    with outfile3 as file:
        for x in zip(*c):
            file.write("{0}\t{1}\n".format(*x))
    outfile3.close()
        
if __name__ == "__main__":
    main()