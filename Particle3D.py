"""
Particle3D; a class to describe point particles in 3D space

An instance describes a particle in Euclidean 3D space: 
(velocity and position are [3] arrays)

Kinetic energy and momentum ([3] array) are calculated for each particle
 
First and second order position updates and velocity update are calculated 
given timestep (and force)

The total kinetic energy, total mass and centre of mass velocity are 
calculated for a system of particles

The positions and velocites for a system of particles are updated. 
 
Separations and closets image separation between two particles are calculated.  

Author: Isaac Watson
 
"""

import numpy as np
import pbc

class Particle3D(object):
    """
    Class to describe point-particles in 3D space

        Properties:
    label: name of the particle
    mass: mass of the particle
    x: x position of the particle
    y: y position of the particle
    z: z position of the particle
    vx: x velocity of the particle
    vy: y velocity of the particle
    vz: z velocity of the particle

        Methods:
    __init__
    __str__
    kinetic_e  - computes the kinetic energy
    momentum - computes the linear momentum
    update_pos_1st - updates the position to 1st order
    update_pos_2nd - updates the position to 2nd order
    update_vel - updates the velocity

        Static Methods:
    new_particle - initializes a P3D instance from a file handle
    sys_kinetic - computes total K.E. of a p3d list
    com_velocity - computes total mass and CoM velocity of a p3d list
    
    """

    def __init__(self, label, mass, x, y, z, vx, vy, vz):
        """
        Initialises a particle in 3D space

        :param label: string, the name of the particle
        :param mass: float, mass of the particle
        :param x: float, x position of particle
        :param y: float, y position of particle
        :param z: float, z position of particle
        :param vx: float, x velocity of particle
        :param vy: float, y velocity of particle
        :param vz: float, z velocity of particle
        
        """
        #define self.parameter so particle properties can be retrieved easily 
        self.name = str(label)
        self.m = float(mass)
        self.p = np.array([x, y, z], float)
        self.v = np.array([vx, vy, vz], float)


    @staticmethod
    def new_particle(i):
        """
        Initialises a Particle3D instance given an input file handle.
            
        :param i: number of particles

        :return Particle3D instance
        
        """
        # Assign particle properties
        label = i
        mass = 1
        x = 0
        y = 0
        z = 0
        vx = 0
        vy = 0
        vz = 0
        
        # Return the particle instance     
        return Particle3D(label, mass, x, y, z, vx, vy, vz)


    def __str__(self):
        # Write position string in the format <label> <x> <y> <z>
        xyz_string = f"{self.name} {self.p[0]} {self.p[1]} {self.p[2]}"
        
        return xyz_string


    def kinetic_e(self):
        """
        Returns the kinetic energy of a Particle3D instance

        :return ke: float, 1/2 m v**2
        
        """
        # Use kinetic energy equation and particle created to caluclate KE as 
        # a float - np.linalg.norm(self.v)=sqrt(vx^2+vy^2+vz^2)
        ke = float(0.5*self.m*(np.linalg.norm(self.v)**2))
        
        return ke


    def momentum(self):
        """
        Returns the linear momentum of a Particle3D instance
        :return p: (3) float np.array, m*v
        
        """
        # Use momentum equation and particle created to calculate p as a float
        p = np.array(self.m*self.v, float)
        
        return p


    def update_pos_1st(self, dt):
        """
        1st order position update

        :param dt: timestep
        
        """
        # Reassign position a new position after timestep
        self.p = self.p+(dt*self.v)


    def update_pos_2nd(self, dt, force):
        """
        2nd order position update

        :param dt: timestep
        :param force: [3] float array, the total force acting on the particle
        
        """
        # Reassign position a new position after timestep
        self.p = self.p+(dt*self.v)+(dt**2*force/(2*self.m))
        

    def update_vel(self, dt, force):
        """
        Velocity update

        :param dt: timestep
        :param force: [3] float array, the total force acting on the particle
        
        """
        # Reassign velocity a new velocity after timestep
        self.v = self.v+(dt*force/self.m)


    @staticmethod
    def sys_kinetic(p3d_list):
        """
        Returns the kinetic energy of the whole system

        :param p3d_list: list in which each item is a P3D instance

        :return sys_ke: \sum 1/2 m_i v_i^2 
        
        """
        # Assign system KE a value of 0 to begin
        sys_ke = 0
        # Loop through each particle in the list
        for i in range(len(p3d_list)):
            # Retrieve KE for each partilce i in the list
            i_ke = p3d_list[i].kinetic_e()
            # Sum particle KEs 
            sys_ke += i_ke

        return sys_ke


    @staticmethod
    def com_velocity(p3d_list):
        """
        Computes the total mass and CoM velocity of a list of P3D's

        :param p3d_list: list in which each item is a P3D instance
        :return total_mass: The total mass of the system 
        :return com_vel: Centre-of-mass velocity
        
        """
        # Assign total mass and centre of mass velocity 0 values to begin 
        total_mass = 0
        com_vel = 0
        # Loop through each particle in the list
        for i in range(len(p3d_list)):
            # Retrive mass of each particle i in the list
            i_mass = p3d_list[i].m
            # Sum retrieved masses
            total_mass += i_mass
            # Retrive momentum of each partilce i in the list
            i_mom = p3d_list[i].momentum()
            # Sum retrieved momentum and divde by total mass
            com_vel += i_mom/total_mass
        
        return total_mass, com_vel
    
    
    @staticmethod
    def update_pos_all(p3d_list, dt, forces_array, l):
        """
        Updates all positions from a list of particles

        :param p3d_list: list in which each item is a P3D instance
        :param dt: timestep of integration loop
        :param force_array: array of forces between particles at that timestep
        :param l: lattice size (boundary condition)
        
        """
        # Loop through particle list updating each position using the 
        # corresponding forces
        for i, p in enumerate(p3d_list):            
            p.update_pos_2nd(dt, forces_array[i])
            p.p = pbc.pbc(p.p, l)
        
    @staticmethod    
    def update_vel_all(p3d_list, dt, forces_array):
        """
        Updates all velocities from a list of particles

        :param p3d_list: list in which each item is a P3D instance
        :param dt: timestep of integration loop
        :param force_array: array of forces between particles at that timestep
        
        """
        # Loop through particle list updating each velocity using the 
        # corresponding forces
        for i, p in enumerate(p3d_list):            
            p.update_vel(dt, forces_array[i])
    
    @staticmethod
    def sep(p1, p2):
        """
        Computes sepration between two particle positions

        :param p1: position of particle 1
        :param p2: position of particle 2
        :return p1-p2: returns particle separation
        
        """
        # Calculate separation and return 
        return p1-p2
    
    @staticmethod
    def mic_sep(p1, p2, l):
        """
        Computes closest image sepration between two particle positions

        :param p1: position of particle 1
        :param p2: position of particle 2
        :param l: lattice size (boundary condition)
        :return mic_sep: returns particle closest image separation
        
        """
        # Calculate the separation of particles in the pbc's
        r = Particle3D.sep(p1, p2)
        r_pbc = pbc.pbc(r, l)
        # Calculate the closest image
        mic_sep = pbc.mic(r_pbc, l)
        
        return mic_sep
    

        
              
        
