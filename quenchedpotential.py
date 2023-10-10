import random
import numpy as np
from wilsongauge import *

def random_su2_R():
    R = np.eye(3)
    R[0][0], R[0][1], R[1][0], R[1][1] = np.random.uniform(-1/2, 1/2, 4)
    return R
def random_su2_S():
    S = np.eye(3)
    S[0][0], S[0][2], S[2][0], S[2][2] = np.random.uniform(-1/2, 1/2, 4)
    return S
def random_su2_T():
    T = np.eye(3)
    T[1][1], T[1][2], T[2][1], T[2][2] = np.random.uniform(-1/2, 1/2, 4)
    return T

N = 10
Nx = Ny = Nt = N

gauge_field_lattice = lattice(Nx, Ny, Nt)

R = 1
T = 3
Nconfig = 1

edges_list = list(gauge_field_lattice.links.keys())


for config in range(Nconfig):
    edge = np.random.choice(edges_list, 1)[0]
    U = gauge_field_lattice.links[edge][3]
    X = random_su2_R() * random_su2_S() * random_su2_T()
    U_candidate = X * U
    print(X)
    '''mc_site = tuple([np.random.randint(0, N), np.random.randint(0, N), np.random.randint(0, N)])
    mc_direction = np.random.randint(1, 3)

    if mc_site[mc_direction] + R < N and mc_site[0] + T < Nt:
        Oi = gauge_field_lattice.wilson_loop(n = mc_site, R = R, j = mc_direction, t = T)

    else:
        continue'''

    

