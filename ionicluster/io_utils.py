import subprocess
from os.path import abspath as path 
from os.path import isdir
from os.path import join
from os import listdir
import numpy as np
import re

class IOUtils:

	"""
	A utility class for file handeling.

	Each of the main methods reports withether it is for InfRETIS, OPES(biased sim), or both.

	
	Attributes
	----------
	current_dir_path : TYPE
	    Description
	load_path : str
	    For InfRETIS, path to the `load` dir, for OPES, path to the dir with trajectories 
	trj_frame_arr : TYPE
	    internal parameter, calculated for infRETIS if sort_to_path is needed
	"""
	
	def __init__(self, load_path : str = 'load'):
		"""
		Initiate the class
		
		Parameters
		----------
		load_path : str, optional
		    path to the load directory
		"""

		self.load_path = load_path


	# ==========================================
	# MAIN METHODS (Public Interface)
	# ==========================================


	def list_paths_dirs(self) -> list:
		"""
		infRETIS
		Finds all the dirs (infRETIS paths) in the "/load" dir.
		
		Returns
		-------
		list
		    the list of directories inside the "load"
		"""
		dir_list = [i for i in listdir(self.load_path) if isdir(join(self.load_path, i)) ]
		dir_list.sort(key=int)

		return dir_list


	def list_trjs(self, current_dir : str) -> list:
		"""
		infRETIS
		Finds all trajectories in the current path (`load/current_dir`).
		Sets the inner parameter `current_dir_path`.
		
		Parameters
		----------
		current_dir : str
		    the name (number) of the current path/directory
		
		Returns
		-------
		list
		    the list of the trajectories
		"""
		self.current_dir_path = join(self.load_path, current_dir, 'accepted')
		return [i for i in listdir(self.current_dir_path) if 'lammpstrj' in i ]


	def find_valid_frames(self, current_dir : str) -> tuple[tuple, int]:
		"""
		infRETIS
		Find all the frames of the current path that were saved in the `traj.txt`.
		Sets the inner parameter `trj_frame_arr`. It is needed for the `sort_to_path` method
		
		Parameters
		----------
		current_dir : str
		    the name (number) of the current path/directory
		
		Returns
		-------
		tuple[tuple, int]
		    valid_frames
		    	-- set of the trj frames that were saved in the `traj.txt`.
		    path_len
		    	-- total lenth of the path
		"""
		trj_dat = np.loadtxt(join(self.load_path, current_dir, 'traj.txt'), dtype=np.dtype('U100'))
		self.trj_frame_arr = trj_dat[:, 1:3]
		valid_frames = set(map(tuple, self.trj_frame_arr)) 

		path_len = len(trj_dat)

		return valid_frames, path_len


	def sort_to_path(self, arr : np.ndarray) -> np.ndarray:
		"""
		infRETIS
		Return the array sorted to reproduce the order of the path from the traj.txt
		
		Parameters
		----------
		arr : np.ndarray
		    the array of interest
		    anything calculated per path step
		    first two colums must contain the trj frame and name like in the traj.txt
		
		Returns
		-------
		np.ndarray
		    sorted array
		"""
		arr_keys = [tuple(row) for row in arr[:, :2]]
		ref_keys = [tuple(row) for row in self.trj_frame_arr]

		indx_order = [arr_keys.index(key) for key in ref_keys]
		return arr[indx_order]


	def list_biased_trjs(self) -> list:
		"""
		OPES/biased simulation
		Finds the list of the OPES trajectories sorted in the chronological order 
		`dump_p.N.lammpstrj` -- this is the expected name of the trj where N identifies
		the chronological sequence 
		
		Returns
		-------
		list
		    the list of the trajectories 
		"""
		def sort_key(s: str) -> list:
		    return [int(p) if p.isdigit() else p for p in re.findall(r'\D+|\d+', s)]

		trjs_list = [i for i in listdir(self.load_path) if 'lammpstrj' in i and 'unwr' not in i ]
		trjs_list.sort(key=sort_key)

		return trjs_list


	def cut_water_from_trjs(self, current_dir : str, N_AT : int):
		"""
		InfRETIS and OPES/biased sim.
		Cut water moleciles from the trjs found with method .list_trjs(current_dir)
		
		Parameters
		----------
		current_dir : str
		    Description
		N_AT : int
		    Total number of CaCO3 atoms 
		    (number of atom after the removal of water molecules per frame)
		"""
		trjs_list = self.list_trjs(current_dir) 

		for trj in trjs_list:

			for trj in trjs_list:
				if 'new' in trj:
					trjs_list.remove(trj)
					subprocess.run(['rm', join(self.current_dir_path, trj)])

				self._cut_water_from_trj(trj, N_AT) #call a helper function, which process a specific trj


	# ==========================================
	# HELPER METHODS (Internal)
	# ==========================================


	def _cut_water_from_trj(self, trj, N_AT):

		new_trj = open(join(self.current_dir_path, 'new_'+trj), 'w')

		with open(join(self.current_dir_path, trj), 'r') as f:

			line = f.readline()

			while line:	


				new_trj.write(line)
				[new_trj.write(f.readline()) for i in range(2)]

				n_at_old = int(f.readline())
				new_trj.write(str(N_AT)+'\n')

				[new_trj.write(f.readline()) for i in range(5)]

				for i in range(n_at_old):

					line = f.readline()
					if len(line.split())<8:
						bad_path.append(d)
						break
					type = line.split()[1]

					if type == '5' or type == '4':
						continue

					else:
						new_trj.write(line)

				line = f.readline()

		new_trj.close()
		subprocess.run(['rm', join(self.current_dir_path, trj)])
		subprocess.run(['mv', join(self.current_dir_path, 'new_'+trj), join(self.current_dir_path, trj)])







		
