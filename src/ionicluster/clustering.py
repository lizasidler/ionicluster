import numpy as np
import MDAnalysis.analysis.distances as distances
from scipy.sparse.csgraph import connected_components
from .structural import StructureAnalyser


class ClusterAnalyser:

	"""
	A class to perform cluster analysys,
	finds atoms belonging to the largest cluster (LC), the CI-cluster, and in the cluster of the specified id.
	Additionally, finds the indices and positions of CI-atoms types A and B separately if pos_B is provided:
		`CI_cl_at_A_ids`, `CI_cl_at_A_pos`, `CI_cl_at_B_ids`, `CI_cl_at_B_pos`. Stores them as attributes.

	This script was employed for the analysis reported in '', where more details on calcluation of each
	quantity is presented.
	
	Attributes
	----------
	box : np.array
	    Simulation box.
	CI_cl_at_A_ids : np.array
	    Indices of CI-cluster atoms of type A. Calculate with the method `CI_cluster` if pos_B is provided.
	CI_cl_at_A_pos : TYPE
	    Positions of CI-cluster atoms of type A. Calculate with the method `CI_cluster` if pos_B is provided.
	CI_cl_at_B_ids : TYPE
	    Indices of CI-cluster atoms of type B. Calculate with the method `CI_cluster` if pos_B is provided.
	CI_cl_at_B_pos : TYPE
	    Positions of CI-cluster atoms of type B. Calculate with the method `CI_cluster` if pos_B is proveded.
	CI_cl_at_ids : TYPE
	    Indices of CI-cluster atoms. Calculate with the method `CI_cluster`.
	CI_cl_at_pos : TYPE
	    Positions of CI-cluster atoms. Calculate with the method `CI_cluster`.
	CI_cl_size : TYPE
	    Number of atoms of types A (and B) in the CI-cluster.
	cl_ids : TYPE
	    Cluster indicies per atom.
	Cluster_by_id_return : TYPE
	    See the description of the method `cluster_by_id`.
	LC_at_A_ids : TYPE
	    Indices of LC atoms of type A. Calculate with the method `largest_cluster` if pos_B is provided.
	LC_at_A_pos : TYPE
	    Positions of LC atoms of type A. Calculate with the method `largest_cluster` if pos_B is provided.
	LC_at_B_ids : TYPE
	    Indices of LC atoms of type B. Calculate with the method `largest_cluster` if pos_B is provided.
	LC_at_B_pos : TYPE
	    Positions of LC atoms of type B. Calculate with the method `largest_cluster` if pos_B is provided.
	LC_at_ids : TYPE
	    Indices of LC atoms.
	LC_at_pos : TYPE
	    Positions of LC atoms.
	LC_size : TYPE
	    Number of atoms of types A (and B) in the LC.
	n_cl : TYPE
	    Number of clusters.
	pos : TYPE
	    Positions of all atoms of interest.
	pos_A : TYPE
	    Positions of all atoms of type A.
	pos_B : TYPE
	    Positions of all atoms of type B.
	unique_cl_ids : TYPE
	    Array if unique cluster ids.
	"""
	
	def __init__(
		self,
		box : np.ndarray| None = None,
		pos_A : np.ndarray| None = None,
		pos_B : np.ndarray | None = None,
		init: bool = False
		):
		"""
		Parameters
		----------
		box : np.ndarray of shape (6,)
		    The dimensions of the simulation box.
        	A 6-element array [lx, ly, lz, alpha, beta, gamma].
		pos_A : np.ndarray of shape (N, 3) | None, optional
		    Coordinates of the first atom group (e.g., central atoms, specific ions, or residues).
		pos_B : np.ndarray of shape (N, 3) | None, optional
		    Coordinates of the secon atom group (e.g., central atoms, specific ions, or residues);
			If None, calculations are done with one atom group only (pos_A against itself).
		init : bool, optional
		    If True, initializes an empty instance of the class
		    without atomic positions or box dimensions. (default is False)
		"""

		if init:
			box = np.array([1, 1, 1, 90, 90, 90]) 
			pos_A = np.array([1, 1, 1]) 

		self.update(box, pos_A, pos_B)


	# ==========================================
	# MAIN METHODS (Public Interface)
	# ==========================================


	def update(
		self,
		box : np.ndarray,
		pos_A : np.ndarray,
		pos_B : np.ndarray | None = None
		):
		"""
		Update the instance with new atomic positions and cell dimensions.
		
		Parameters
		----------
		box : np.ndarray of shape (6,)
		    The dimensions of the simulation box.
        	A 6-element array [lx, ly, lz, alpha, beta, gamma].
		pos_A : np.ndarray of shape (N, 3) | None, optional
		    Coordinates of the first atom group (e.g., central atoms, specific ions, or residues).
		pos_B : np.ndarray of shape (N, 3) | None, optional
		    Coordinates of the secon atom group (e.g., central atoms, specific ions, or residues);
			If None, calculations are done with one atom group only (pos_A against itself).
		
		Raises
		------
		ValueError
		    Expected array shape (6,), but got {box.shape}
		"""
		if box.shape != (6,):
		    raise ValueError(f"Expected array shape (6,), but got {box.shape}")

		self.box = box
		self.pos_A = pos_A
		self.pos_B = pos_B

		if pos_B is not None:
		    # If the atom counts don't change frame-to-frame, 
		    # using slicing assignment here prevents re-allocation!
		    if hasattr(self, 'pos') and self.pos.shape == (pos_A.shape[0] + pos_B.shape[0], 3):
		        self.pos[:pos_A.shape[0]] = pos_A
		        self.pos[pos_A.shape[0]:] = pos_B
		    else:
		        # Runs the first time (or if array sizes dynamically change)
		        self.pos = np.vstack([pos_A, pos_B])
		else:
		    self.pos = pos_A

		self.cl_ids, self.n_cl = None, None


	def count(self, rcut : float = 4.0) -> tuple[int, np.ndarray]:
		"""
		Find clusters of atoms within a specified cutoff distance.

		Update the internal array of cluster IDs per atom (`cl_ids`), unique cluster 
		IDs (`unique_cl_ids`), and the total number of clusters (`n_cl`).
		
		Parameters
		----------
		rcut : float, optional
		    cutoff distance
		
		Returns
		-------
		tuple[np.ndarray, np.ndarray]
		    n_cl -- Number of clusters
		    cl_ids -- The cluster labels per atom.
		"""

		cont_matr = distances.contact_matrix(self.pos, cutoff=rcut, returntype='sparse', box=self.box) 
		self.n_cl, self.cl_ids = connected_components(csgraph=cont_matr, directed=False, return_labels=True)

		self.unique_cl_ids = set(self.cl_ids)

		return self.n_cl, self.cl_ids


	def CI_cluster(
		self,
		cl_ids : np.ndarray | None = None,
		central_at_id : int | None = None
		) -> tuple[np.ndarray, np.ndarray, int]:
		"""
		Calculate the indices and positions of atoms belonging to the CI-cluster (cluster containing cental ion):
		(`CI_cl_at_ids`), (`CI_cl_at_pos`) 
		and the CI-cluster size (`CI_cl_size`).

		Additionally, finds the indices and positions of CI-atoms types A and B separately if pos_B is provided:
		`CI_cl_at_A_ids`, `CI_cl_at_A_pos`, `CI_cl_at_B_ids`, `CI_cl_at_B_pos`. Stores them as attributes.
		
		The ID of the central atom will be calculated if not provided.
		
		Parameters
		----------
		cl_ids : np.ndarray | None, optional
		    The cluster labels per atom (e.g., as returned by the `.count()` method) may be provided directly.
	    	If None, uses the internal `cl_ids` attribute or evaluates them in-place.
		central_at_id : int | None, optional
		    The ID of the central atom relative to the `pos_A` array. If None, 
		    the central atom will be determined internally using method .x_min_distance() 
			from Structural Analysis.
		
		Returns
		-------
		tuple[np.ndarray, np.ndarray, int]
		    -- np.ndarray
		    An array containing the positions of the atoms in the CI-cluster.
		    -- np.ndarray
		    An array containing the IDs of the atoms in the CI-cluster.
		    -- int
		    The total number of atoms belonging to the CI-cluster.
		"""

		self._set_cl_ids(cl_ids) # helper function that organises cl_ids

		if central_at_id is None:

			struct = StructureAnalyser(self.box, pos_A=self.pos_A, pos_B=self.pos_B)
			_, central_at_id = struct.x_min_distance()

		CI_cl_id = self.cl_ids[central_at_id]
		self.CI_cl_at_ids = np.where(self.cl_ids==CI_cl_id)[0]
		self.CI_cl_size = len(self.CI_cl_at_ids)
		self.CI_cl_at_pos = self.pos[self.CI_cl_at_ids]

		if self.pos_B is not None:
			self.CI_cl_at_A_ids = self.CI_cl_at_ids[self.CI_cl_at_ids<len(self.pos_A)]
			self.CI_cl_at_A_pos = self.pos_A[self.CI_cl_at_A_ids]
			self.CI_cl_at_B_ids = self.CI_cl_at_ids[self.CI_cl_at_ids>=len(self.pos_A)] - len(self.pos_A)
			self.CI_cl_at_B_pos = self.pos_B[self.CI_cl_at_B_ids]

		return self.CI_cl_at_pos, self.CI_cl_at_ids, self.CI_cl_size


	def largest_cluster(self, cl_ids : np.ndarray | None = None
		) -> tuple[np.ndarray, np.ndarray, int]:
		"""
		Calculate the indices and positions of atoms belonging to the LC (largest cluster):
		(`LC_at_ids`), (`LC_at_pos`)
		and the LC size (`LC_size`).

		Additionally, finds the indices and positions of LC atoms types A and B separately if pos_B is provided:
		`LC_at_A_ids`, `LC_at_A_pos`, `LC_at_B_ids`, `LC_at_B_pos`. Stores them as attributes.

		
		Parameters
		----------
		cl_ids : np.ndarray | None, optional
		    The cluster labels per atom (e.g., as returned by the `.count()` method) may be provided directly.
	    	If None, uses the internal `cl_ids` attribute or evaluates them in-place.
		
		Returns
		-------
		tuple[np.ndarray, np.ndarray, int]
			-- np.ndarray
		    An array containing the positions of the atoms in the LC.
		    -- np.ndarray
		    An array containing the IDs of the atoms in the LC.
		    -- int
		    The total number of atoms belonging to the LC.
		"""
		self._set_cl_ids(cl_ids) # helper function that organises cl_ids

		self.LC_at_ids = np.where(self.cl_ids==np.argmax(np.bincount(self.cl_ids)))[0]
		self.LC_size = len(self.LC_at_ids)
		self.LC_at_pos = self.pos[self.LC_at_ids]

		if self.pos_B is not None:
			self.LC_at_A_ids = self.LC_at_ids[self.LC_at_ids<len(self.pos_A)]
			self.LC_at_A_pos = self.pos_A[self.LC_at_A_ids]
			self.LC_at_B_ids = self.LC_at_ids[self.LC_at_ids>=len(self.pos_A)] - len(self.pos_A)
			self.LC_at_B_pos = self.pos_B[self.LC_at_B_ids]

		return self.LC_at_pos, self.LC_at_ids, self.LC_size

	Cluster_by_id_return = (
	    tuple[np.ndarray, np.ndarray, int] | 
	    tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]
	)


	def cluster_by_id(
		self,
		cl_id,
		cl_ids : np.ndarray | None = None,
		return_AB_sep : bool = False
		) -> Cluster_by_id_return:
		"""
		Calculate the indices and positions of atoms belonging to the cluster with the given `cl_id`:
		(`cl_at_ids`), (`cl_at_pos`).
		and the size of this cluster (`cl_size`).

		If return_AB_sep, find and returns the indices and positions of cluster aroms types A and B separately
		if pos_B is provided:
		`cl_at_A_ids`, `cl_at_A_pos`, `cl_at_B_ids`, `cl_at_B_pos`.
		
		Parameters
		----------
		cl_id : np.ndarray | None, optional
		    The cluster labels per atom (e.g., as returned by the `.count()` method) may be provided directly.
	    	If None, uses the internal `cl_ids` attribute or evaluates them in-place.
	    cl_ids : np.ndarray | None, optional
		    The cluster labels per atom (e.g., as returned by the `.count()` method) may be provided directly.
	    	If None, uses the internal `cl_ids` attribute or evaluates them in-place.
	    return_AB_sep : bool, optional
	    	Determins whether to calculate the indices and positions of cluster aroms types A and B separately
		
		Returns
		-------
		Cluster_by_id_return:
			if not return_AB_sep:
			    tuple[np.ndarray, np.ndarray, int]
			    -- np.ndarray
			    An array containing the pos of the cl_id cluster atoms.
			    -- np.ndarray
			    An array containing the IDs of the atoms in the cl_id cluster.
			    -- int
			    The total number of atoms belonging to the cl_id cluster.

			else:
				tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]
			    -- np.ndarray
			    An array containing the pos of the cl_id cluster atoms.
			    -- np.ndarray
			    An array containing the IDs of the atoms in the cl_id cluster.
			    -- np.ndarray
			    An array containing the pos of the cl_id cluster atoms type A.
			    -- np.ndarray
			    An array containing the IDs of the type A atoms in the cl_id cluster.
			    -- np.ndarray
			    An array containing the pos of the cl_id cluster atoms type B.
			    -- np.ndarray
			    An array containing the IDs of the type B atoms in the cl_id cluster.
			    -- int
			    The total number of atoms belonging to the cl_id cluster.

		"""
		self._set_cl_ids(cl_ids) # helper function that organises cl_ids

		cl_at_ids=np.where(self.cl_ids==cl_id)[0]
		cl_size = len(cl_at_ids)
		cl_at_pos = self.pos[cl_at_ids]

		if self.pos_B is not None and return_AB_sep:
			cl_at_A_ids = cl_at_ids[cl_at_ids<len(self.pos_A)]
			cl_at_A_pos = self.pos_A[cl_at_A_ids]
			cl_at_B_ids = cl_at_ids[cl_at_ids>=len(self.pos_A)] - len(self.pos_A)
			cl_at_B_pos = self.pos_B[cl_at_B_ids]

			return cl_at_pos, cl_at_ids, cl_at_A_pos, cl_at_A_ids, cl_at_B_pos, cl_at_B_ids, cl_size

		return cl_at_pos, cl_at_ids, cl_size


	@staticmethod
	def find_ox_ids(cl_at_c_id : np.array) -> np.array:
		"""
		Calculate oxygen atom ids (Oxygen from CO3 ion) if ids of carbon atoms are known.
		
		Parameters
		----------
		cl_at_c_id : np.array
		    Ids of carbon atoms (part of CO3 ions).
		
		Returns
		-------
		np.array
		    Oxygen atom ids.
		"""
		cl_at_o_id = np.zeros(len(cl_at_c_id)*3, dtype=np.int64)

		for i,_id in enumerate(cl_at_c_id):
			i0 = i*3
			_id0 = _id*3

			cl_at_o_id[i0] = _id0
			cl_at_o_id[i0+1] = _id0+1
			cl_at_o_id[i0+2] = _id0+2


		return cl_at_o_id


	# ==========================================
	# HELPER METHODS (Internal)
	# ==========================================


	def _set_cl_ids(self, cl_ids : np.ndarray | None = None):
		"""
		Helper function
		If cl_ids is provided,
		Update class attributes `cl_ids`, `unique_cl_ids`, and `n_cl`.

		Otherwise, evaluates these attributes with .count()
		
		Parameters
		----------
		cl_ids : np.array | None, optional
		    The cluster labels per atom.
		"""
		if all(v is None for v in[cl_ids, self.cl_ids, self.n_cl]):

			self.n_cl, self.cl_ids = self.count(rcut = 4.0)

		elif cl_ids is not None:

			self.cl_ids = cl_ids
			self.unique_cl_ids = set(cl_ids)
			self.n_cl = len(self.unique_cl_ids)



