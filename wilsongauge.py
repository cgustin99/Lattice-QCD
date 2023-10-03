import numpy as np
from itertools import product


#START IN 3D

#Paramters
Nx, Ny, Nt = 3, 3, 1
N_links_mu0 = Nx * Ny * Nt
N_links_mu1 = (Nx - 1) * Ny * Nt
N_links_mu2 = Nx * (Ny - 1) * Nt



class qcd_lattice():
    def __init__(self, Nx, Ny, Nt): 
        self.Nx = Nx
        self.Ny = Ny
        self.Nt = Nt

        lattice = np.zeros((Nx, Ny, Nt))

Umu0 = np.kron(np.eye(N_links_mu0), np.eye(3)) #Initial lattice configuration 
Umu1 = np.kron(np.eye(N_links_mu1), np.eye(3))
Umu2 = np.kron(np.eye(N_links_mu2), np.eye(3))
Umu = np.array([Umu0, Umu1, Umu2])

full_lattice = qcd_lattice(Nx, Ny, Nt)

def wilson_gauge_action(lattice, U):
    #Compute S_G(U) = 2/g^2 \sum_{n \in \Lambda} \sum_{\mu < \nu} Re tr[1 - U_mu,nu(n)]
    S_G = 0

    for n in product(range(lattice.Nx), range(lattice.Ny), range(lattice.Nt)):

        for mu in (0, 1, 2):
            for nu in range(mu, 2):
                Uuv = 1#Plaquette
                
    return

wilson_gauge_action(full_lattice, Umu0)