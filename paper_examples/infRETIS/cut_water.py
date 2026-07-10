from ionicluster import IOUtils
import MDAnalysis as mda

"""

Removes water molecules from the path trajectories (to save disk space)

"""

u = mda.Universe('chi_0.01_nion_72_cac_OP/lammps_input/lammps.data', format='DATA')

Ca = u.select_atoms('type 1')
C = u.select_atoms('type 2')
O = u.select_atoms('type 3')
N_CaCO3 = len(Ca) + len(C) + len(O)

util = IOUtils('chi_0.01_nion_72_cac_OP/load')
dirs = util.list_paths_dirs()

for current_dir in dirs:

	util.cut_water_from_trjs(current_dir, N_CaCO3)