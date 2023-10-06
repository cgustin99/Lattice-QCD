import numpy as np
from wilsongauge import *

Nx, Ny, Nt = 10, 10, 10

lattice = qcd_lattice(Nx, Ny, Nt)


def mcmc_pathintegral():
    #Approximates <O> = \int D[Ψ, \barΨ, U] O[Ψ, barΨ, U] exp(-S)
    
    return