# -*- coding: utf-8 -*-
"""
A set of functions to caluclate the forc ena dpotneital of a sytem of 
particles in a Lennard Jones potential 
@author: Isaac Watson
"""
import numpy as np
from Particle3D import Particle3D


def force_lj(sep, cutoff):
    """
    Computes the force between 2 particles in a Lennard Jones potential
    
    :param sep : separation between two particles (vector)
    :param cutoff : max separation to complete calculations (float)
    :return force: force between 2 particles (float)

    """
    # Magnitude of the separation
    r = np.linalg.norm(sep)
    # Calculate force if separation is less than cutoff of vector
    if r < cutoff:
        force = 48*(r**-14 - 0.5*r**-8)*sep
    else:
        force = 0
    
    return force

def get_forces(positions, l, cutoff): 
    """
    Creates an array of all forces between particles
    
    :param positions : current positions (array)
    :param l : lattice size (vector)
    :param cutoff : max distance to complete calculations (float)
    :return p_forces : forces on all particles (array)

    """
    N = len(positions)
    # Create empty array of all forces
    p_forces = np.zeros([N, 3])
    # Get force of particle j on i
    for i in range(N):
        for j in range(i):
            if i != j:
                # Vector (3)
                micsep = Particle3D.mic_sep(positions[i], positions[j], l)
                # Vector (3)
                force = force_lj(micsep, cutoff)
                # Array element = vector 
                p_forces[i] += force
                # Force on particle j due to Newtons 3rd Law
                p_forces[j] -= force       
    
    return p_forces


def potential_lj(sep, cutoff): 
    """
    :param sep : separation between two particles (vector)
    :param cutoff : max distance to complete calculations (float)
    :return potential : value for potential between two particles (float)

    """
    # Magnitude of the separation
    r = np.linalg.norm(sep)  
    # Lennard Jones potential if separation is less than cutoff
    if r < cutoff:
        potential = 4*(r**-12-r**-6)-4*(cutoff**-12-cutoff**-6)
    else:
        potential = 0
    
    return potential 

def sys_potential_lj(positions, l, cutoff): 
    """
    Computes the potential of the entire system of particles 
    
    :param positions : current positions (array)
    :param cutoff : max distance to complete calculations (float)
    :return sys_pot : potential of the entire system of particles (float)

    """
    # Number of particles
    N = len(positions)
    # Set intial potential energy for no particles
    sys_pot = 0 
    # Loop through all particles and calculate the potential energy 
    for i in range(N):
        for j in range(i):
            if i != j:
                micsep = Particle3D.mic_sep(positions[i], positions[j], l)
                sys_pot += potential_lj(micsep, cutoff)
                            
    return sys_pot