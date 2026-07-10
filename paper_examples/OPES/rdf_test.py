import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from os.path import join

from ionicluster import StructureAnalyser
from ionicluster import ClusterAnalyser
from ionicluster import IOUtils
"""

Calculate RDF (PDF) of the CI-cluster (Ca central) using the OPES trajectories

"""

bins = 250 

rdf_Ca_C = np.zeros(bins)
rdf_Ca_Ca = np.zeros(bins)
rdf_C_C = np.zeros(bins)
rdf_Ca_O = np.zeros(bins)
rdf_all = np.zeros(bins)

count_ca_c = 0
count_ca_ca = 0
count_c_c = 0
count_ca_o = 0
count_all = 0

# set up analysers

cl = ClusterAnalyser(init = True)
_ca = StructureAnalyser(init = True)
_c = StructureAnalyser(init = True)
_all = StructureAnalyser(init = True)

io_utils = IOUtils('.')

# if an OPES trajectory is split across multiple files via restarts 
# find all OPES trajectories in chronological order

trjs_list = io_utils.list_biased_trjs()


for trj in trjs_list[:1]:

	u = mda.Universe(join(io_utils.load_path, trj), format='LAMMPSDUMP')

	Ca = u.select_atoms('type 1')
	C = u.select_atoms('type 2')
	O = u.select_atoms('type 3')

	for ts in u.trajectory[:99]:

			box = u.dimensions
			pos_A = Ca.positions
			pos_B = C.positions
			pos_O = O.positions

			cl.update(box, pos_A, pos_B)

			# cluster analysis for the current frame; find CI-cluster

			cl.CI_cluster()

			CI_pos_ca = cl.CI_cl_at_A_pos
			CI_pos_c = cl.CI_cl_at_B_pos

			o_cl_id_arr = cl.find_ox_ids(cl.CI_cl_at_B_ids)
			CI_pos_o = pos_O[o_cl_id_arr] # change in the case of the C-Ca 6th min OP 

			_ca.update(box, CI_pos_ca)
			_c.update(box, CI_pos_c)
			_all.update(box, np.vstack([CI_pos_ca, CI_pos_c, CI_pos_o]))

			# calculate RDF for the CI-cluster
			
			if CI_pos_ca.shape[0] and CI_pos_c.shape[0]:
		
				rdf, bin_centers = _ca.rdf(pos_B = CI_pos_c)
				rdf_Ca_C += rdf
				count_ca_c += 1

				rdf, bin_centers = _all.rdf()
				rdf_all += rdf
				count_all += 1

				rdf, bin_centers = _ca.rdf(pos_B = CI_pos_o)
				rdf_Ca_O += rdf
				count_ca_o += 1


			if CI_pos_ca.shape[0] > 1:
			    rdf, bin_centers = _ca.rdf()
			    rdf_Ca_Ca += rdf
			    count_ca_ca += 1

			if CI_pos_c.shape[0] > 1:
				rdf, bin_centers = _c.rdf()
				rdf_C_C += rdf
				count_c_c += 1
	


rdf_avg_Ca_C = rdf_Ca_C / count_ca_c
rdf_avg_Ca_O = rdf_Ca_O / count_ca_o
rdf_avg_Ca_Ca = rdf_Ca_Ca / count_ca_ca
rdf_avg_C_C = rdf_C_C / count_c_c
rdf_avg_all = rdf_all / count_all

plt.figure(figsize=(10, 6))

plt.plot(bin_centers, rdf_avg_Ca_C/rdf_avg_Ca_C[-1], label='Ca-C', linewidth=2)
plt.plot(bin_centers, rdf_avg_Ca_O/rdf_avg_Ca_O[-1], label='Ca-O', linewidth=2)
plt.plot(bin_centers, rdf_avg_Ca_Ca/rdf_avg_Ca_Ca[-1], label='Ca-Ca', linewidth=2)
plt.plot(bin_centers, rdf_avg_C_C/rdf_avg_C_C[-1], label='C-C', linewidth=2)
plt.plot(bin_centers, rdf_avg_all/rdf_avg_all[-1], label='all', linewidth=2)

plt.xlabel('Distance (Å)', fontsize=12)
plt.ylabel('g(r)', fontsize=12)
plt.title('Radial Distribution Functions', fontsize=14)
plt.legend(fontsize=11)
plt.grid(True, alpha=0.3)
plt.tight_layout()

#plt.savefig('rdf_plots_all.png', dpi=300, bbox_inches='tight')
plt.show()
