from pydoc import locate
import numpy as np
from itertools import product
from collections import defaultdict


#START IN 3D

#Paramters
Nx, Ny, Nt = 3, 3, 3
g = 1 #Coupling
nDim = 3 #Number of dimensions

def I3():
    return np.eye(3)

def euclidean_distance(n1, n2):
    return np.sqrt((n1[0] - n2[0])**2 + (n1[1] - n2[1])**2 + (n1[2] - n2[2])**2)

def locate_link_dir(ni, nj):
    minus = [np.abs(a_i - b_i) for a_i, b_i in zip(ni, nj)]
    return minus.index(1)

def test_valid_pt(n):
    if n[0] > Nt - 1:
        n[0] -= 2
    if n[1] > Nx - 1:
        n[1] -= 2
    if n[2] > Ny - 1:
        n[2] -= 2
    return n

def wilson_line_path(n, R, j):
        path = [n]

        for k in range(1, R + 1 - n[j]):
            var_link = list(path[-1])
            var_link[j] += 1
            path.append(tuple(var_link))

        return path

def wilson_line_path_dagger(path, R, j):
    path_dag = []
    for site in path:
        var_site = list(site)
        var_site[j] += R
        path_dag.append(tuple(var_site))
    return path_dag[::-1] 

def tidy_path(path):
    tidied = [path[0]]

    for site in path:
        if site != tidied[-1]:
            tidied.append(site)
    return tidied
    

class qcd_lattice():
    def __init__(self, Nx, Ny, Nt): 

        self.Nx = Nx
        self.Ny = Ny
        self.Nt = Nt

        nlinks = 0
        self.links = {}
        #links = {link  #: [sites associated with this link], direction of link, corresponding U}

        self.states = list(product(range(0, Nt), range(0, Nx), range(0, Ny)))
        
        self.sites = dict((site, []) for site in self.states)
        #sites = {site: {associated link: direction link points}}

        for i in range(len(self.states)):
            for j in range(i, len(self.states)):
                ni, nj = self.states[i], self.states[j]
                if euclidean_distance(ni, nj) == 1.0: 
                    nlinks += 1 
                    direction = locate_link_dir(ni, nj)
                    self.links[nlinks] = [ni, nj, direction, I3()]
                    self.sites[ni].append(nlinks), self.sites[nj].append(nlinks)

        self.Uinit = np.kron(np.eye(nlinks), np.eye(3))

    def get_shared_edge(self, n, m):
        #m and n must be nearest neighbors!
        n_edges, m_edges = self.sites[tuple(n)], self.sites[tuple(m)]
        shared_edge = list(set(n_edges) & set(m_edges))[0]
        return shared_edge

    def Plaquette(self, n, mu, nu):
        plq = np.eye(3)

        n_var = list(n)

        n_plus_mu, n_plus_nu, n_plus_mu_nu = list(n_var), list(n_var), list(n_var)

        n_plus_mu[mu] += 1
        n_plus_nu[nu] += 1
        n_plus_mu_nu[mu] += 1
        n_plus_mu_nu[nu] += 1


        plq_sites = [test_valid_pt(n_var), test_valid_pt(n_plus_mu), 
            test_valid_pt(n_plus_mu_nu), test_valid_pt(n_plus_nu), n_var]

        path = 1
        for ni, nj in zip(plq_sites, plq_sites[1:]):
            shared_edge = self.get_shared_edge(ni, nj)
            if path == 1 or path == 2:
                plq *= self.links[shared_edge][3] 
            elif path == 3 or path == 4:
                plq *= np.conjugate(self.links[shared_edge][3])
            
        return plq

    
    
    def wilson_loop(self, n, R, j, t, temp_gauge = False):
        #n: Initial Point
        #R: spatial distance between m and n 
        #j: Direction
        #t: Temporal distance (from t = 0 to t = t)

        loop = np.eye(3)
        num_edges = 0

        assert R < self.Nx, "Loop goes beyond bounds of the lattice"
        assert t < self.Ny, "Loop goes beyond bounds of the lattice"
        

        S = wilson_line_path(n, R, j)
        T = wilson_line_path(S[-1], t, 0)
        Sdag = wilson_line_path_dagger(S, t, 0)
        Tdag = wilson_line_path_dagger(T, -R, j)

        if temp_gauge == True: T, Tdag = np.eye(3), np.eye(3)
        
        path = tidy_path(S + T + Sdag + Tdag)

        for site_i in range(len(path) - 1):
            shared_edge = self.get_shared_edge(path[site_i], path[site_i + 1])
            num_edges += 1
            if num_edges < R + t:
                loop *= self.links[shared_edge][3]
            else: loop *= np.conjugate(self.links[shared_edge][3])

        
        return np.trace(loop)


full_lattice = qcd_lattice(Nx, Ny, Nt)

def wilson_gauge_action(lattice):
    #Compute S_G(U) = 2/g^2 \sum_{n \in \Lambda} \sum_{\mu < \nu} Re tr[1 - U_mu,nu(n)]
    S_G = 0
    beta = 6/(g**2)

    for n in lattice.states:
        for mu in (0, 1, 2):
            for nu in range(mu + 1, 3):
                Uuv = lattice.Plaquette(n, mu, nu) #Plaquette
                S_G += np.real(np.trace(np.eye(3) - Uuv))
    S_G *= beta

    return S_G

