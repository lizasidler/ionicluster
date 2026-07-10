import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from os.path import join

from ionicluster import StructureAnalyser
from ionicluster import ClusterAnalyser
from ionicluster import IOUtils
"""

Calculate Rg for each cluster 
and combine them into the Rg distribution per cluster size:

rg_hist_per_cl_size, shape (n_caco3+1, bins)

"""

n_caco3 = 72

io_utils = IOUtils('chi_0.01_nion_72_cac_OP/load')
dirs = io_utils.list_paths_dirs()

rg_hist_per_cl_size = np.zeros((n_caco3+1, 250), dtype=float)
rg_array = np.linspace(0, 15, 250)

cl = ClusterAnalyser(init = True)
struct = StructureAnalyser(init = True)

for d in tqdm(dirs):

	valid_frames, _ = io_utils.find_valid_frames(d)
	trjs_list = io_utils.list_trjs(d)
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

			# full cluster analysis

			cl.update(box, ca_pos, c_pos)
			struct.update(box)

			cl.count()

			# go via all clusters

			for cl_id in cl.unique_cl_ids:

				cl_at_pos, _, cl_size = cl.cluster_by_id(cl_id=cl_id)

				if cl_size < 3: 
					continue

				# calculate Rg, for the current cluster

				Rg, _, _, _, _, _ = struct.Rg_K1_K2(pos=cl_at_pos)

				rg_id = np.argmin(abs(rg_array-Rg))
				rg_hist_per_cl_size[cl_size, rg_id]+=1



#np.savetxt('rg_hist_cl_size.txt', rg_hist_per_cl_size, fmt="%s")

plt.figure(figsize=(6,4))
plt.plot(rg_array, rg_hist_per_cl_size[18]) # plor Rg of the clusters with size 18

#splt.savefig('rg_12.png')
plt.show()



