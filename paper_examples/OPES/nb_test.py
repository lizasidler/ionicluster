import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from os.path import join

from ionicluster import StructureAnalyser
from ionicluster import IOUtils

"""

Calculate CN for each atom

"""

io_utils = IOUtils('.')

# if an OPES trajectory is split across multiple files via restarts 
# find all OPES trajectories in chronological order

trjs_list = io_utils.list_biased_trjs()

_ca = StructureAnalyser(init = True)
_c = StructureAnalyser(init = True)

file = open('nb_list.txt', 'w')

for trj in trjs_list:

	u = mda.Universe(join(io_utils.load_path, trj), format='LAMMPSDUMP')

	Ca = u.select_atoms('type 1')
	C = u.select_atoms('type 2')

	for ts in u.trajectory:
	
		box = u.dimensions
		ca_pos = Ca.positions
		c_pos = C.positions

		_ca.update(box, ca_pos)
		_c.update(box, c_pos)

		# calculate different CN per atom

		ca_c_matr = _ca.distance_matrix(pos_B = c_pos)

		nb_ca_c = _ca.CN(dist=ca_c_matr, d0=3.5, dmax=4.5)
		nb_ca_ca = _ca.CN(ca_pos, d0=4.5, dmax=5.5)

		nb_c_ca = _c.CN(dist=ca_c_matr.T, d0=3.5, dmax=4.5)
		nb_c_c = _c.CN(d0=5.1, dmax=6.1)

		# dump CN per atom in the file (one file across all trajectories)

		for i in range(len(ca_pos)):

			file.write(str(nb_ca_c) + ' ' + str(nb_ca_ca) + ' ' + str(nb_c_ca) + ' ' + str(nb_c_c)+' ')
			file.write('\n')

file.close()




	
