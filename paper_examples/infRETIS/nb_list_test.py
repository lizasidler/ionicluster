import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from os.path import join

from ionicluster import StructureAnalyser
from ionicluster import IOUtils

"""

Calculate CN for each atom for each infRETIS path
this script would generate new file in each dir

"""


n_units = int(72/2) #for the current study, number of the CaCO3 units

io_utils = IOUtils('chi_0.01_nion_72_cac_OP/load')
dirs = io_utils.list_paths_dirs()

_ca = StructureAnalyser(init = True)
_c = StructureAnalyser(init = True)

for d in tqdm(dirs):

	valid_frames, path_len = io_utils.find_valid_frames(d)
	trjs_list = io_utils.list_trjs(d)

	# new CN array for each dir
	
	cn_array = np.zeros((path_len, n_units*4+2), dtype=np.dtype('U100'))
	step = 0

	for trj in trjs_list:

		u = mda.Universe(join(io_utils.load_path, d, 'accepted', trj), format='LAMMPSDUMP')

		Ca = u.select_atoms('type 1')
		C = u.select_atoms('type 2')

		for ts in u.trajectory:

			if (trj, str(ts.frame)) not in valid_frames:
				continue  # Skip because it wasn't saved in traj.txt!

			box = u.dimensions
			ca_pos = Ca.positions
			c_pos = C.positions

			cn_array[step, 0] = trj
			cn_array[step, 1] = str(ts.frame)

			_ca.update(box, ca_pos)
			_c.update(box, c_pos)

			# calculate different CN per atom

			ca_c_matr = _ca.distance_matrix(pos_B = c_pos)

			cn_ca_c = _ca.CN(dist=ca_c_matr, d0=3.5, dmax=4.5)
			cn_ca_ca = _ca.CN(ca_pos, d0=4.5, dmax=5.5)

			cn_c_ca = _c.CN(dist=ca_c_matr.T, d0=3.5, dmax=4.5)
			cn_c_c = _c.CN(d0=5.1, dmax=6.1)

			# write all CNs into the array corresponding to the current dir

			for i in range(len(ca_pos)):

				cn_array[step, 2 + i*4: 2 + (i+1)*4] = cn_ca_c[i], cn_ca_ca[i], cn_c_ca[i], cn_c_c[i]

			step += 1

	cn_array = io_utils.sort_to_path(cn_array) # sort array to follow the path (based on traj.txt)
	np.savetxt(join(io_utils.load_path, str(d), 'cn_per_atom.txt'), cn_array, fmt="%s")



	
