import numpy as np
from wilsongauge import *

N = 10
Nx = Ny = Nt = N

gauge_lattice = qcd_lattice(Nx, Ny, Nt)

R = 1

mc_site = [np.random.randint(0, N), np.random.randint(0, N), np.random.randint(0, N)]
mc_direction = np.random.randint(0, 3)


