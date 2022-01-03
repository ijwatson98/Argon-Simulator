"""
Computer Modelling Excercise 1 - PBC and MIC functions

Author: Isaac Watson

"""

import numpy as np

def pbc(pos, l):
    """
    Place particle within pbc's
    
    :param pos : particle position
    :param l : lattice size
    :return :position in pbc

    """
    # Use mod function to find position of particle within pbc
    return np.mod(pos, l)

def mic(pos, l):
    """
    Retrieve closest image of the particle

    :param pos : particle posisition
    :param l : lattice size
    :return : nearest image particle once in pbc

    """
    # Add half the length of the box to each vector component and use mod 
    # function to return new coord to pbc then take back away half the length 
    # of the box to produce coordinates of the closest image
    return np.mod(pos+l/2, l) - l/2


