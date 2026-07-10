import numpy as np
import MDAnalysis.analysis.distances as distances

class StructureAnalyser:


	"""
	A class for geometric and structural analysis of atomic configurations.

    This class provides methods to calculate distance matrices, coordination
    numbers (CN), radii of gyration (Rg), shape parameters (K1, K2),
    radial (partial) distribution finctions (RDF), and
    x minimal distances under periodic boundary conditions.

    This script was employed for the analysis reported in '', where more details on calcluation of each
	quantity is presented.
	
	Attributes
	----------
	box : np.ndarray
	    Description
	pos_A : np.ndarray
	    Coordinates of the first atom group (Group A).
	pos_B : np.ndarray
	    Coordinates of the second atom group (Group B).
	"""
	
	def __init__(
		self,
		box : np.ndarray| None = None,
		pos_A : np.ndarray | None = None,
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
		
		Raises
		------
		ValueError
		    "Expected array shape (6,), but got {box.shape}"
		"""
		if init:
			box = np.array([1, 1, 1, 90, 90, 90])  

		self.update(box, pos_A, pos_B)


	def update(
		self,
		box : np.ndarray,
		pos_A : np.ndarray | None = None,
		pos_B : np.ndarray | None = None,
		):

			if box.shape != (6,):
				raise ValueError(f"Expected array shape (6,), but got {box.shape}")

			self.box = box
			self.pos_A = pos_A
			self.pos_B = pos_B
		
	# ==========================================
	# MAIN METHODS (Public Interface)
	# ==========================================


	def Rg_K1_K2(
		self,
		pos : np.ndarray,
		masses : np.ndarray | None = None
		) -> tuple[float, float, float, float, float, float]:
		"""
		Calculate the radius of gyration (Rg), shape parameters (K1, K2), 
		and eigenvalues of the tensor of gyration (L1, L2, L3).

		Parameters
		----------
		pos : np.ndarray of shape (N, 3)
		    The coordinates of the atoms.
		masses : np.ndarray | None, optional
		    The masses of the atoms used for center-of-mass weighting. 
	        If None, all atoms are assumed to have equal mass, 
	        geometric ragius of gyration is calculated.
		
		Returns
		-------
		tuple[float, float, float, float, float, float]
		    Rg -- radius of gyration.
		    K1, K2 -- shapes parameters
		    L1, L2, L3 -- eigenvalues of the tensor of gyration.

		"""

		if masses is None:
			masses = np.ones(len(pos))

		# calculate center of mass (com) considering PBC

		scaled_pos = pos / self.box[:3]

		u = np.cos(2*np.pi*scaled_pos)
		v = np.sin(2*np.pi*scaled_pos)

		total_mass = np.sum(masses)

		u_com = np.average(u, weights=masses, axis=0)
		v_com = np.average(v, weights=masses, axis=0)

		com = (1 / (2 * np.pi) * np.arctan2(v_com, u_com)) * self.box[:3]

		# calculate tensor of gyration

		dist = self._pbc(pos - com)

		M = sum(masses)

		rg_xx = sum(dist[:, 0]**2*masses)/M
		rg_yy = sum(dist[:, 1]**2*masses)/M
		rg_zz = sum(dist[:, 2]**2*masses)/M

		rg_xy = sum(dist[:, 0]*dist[:, 1]*masses)/M
		rg_xz = sum(dist[:, 0]*dist[:, 2]*masses)/M
		rg_yz = sum(dist[:, 1]*dist[:, 2]*masses)/M

		rg_ten = np.array([ rg_xx, rg_yy, rg_zz, rg_xy, rg_xz, rg_yz ])

		matrix = np.array([[rg_ten[0], rg_ten[3], rg_ten[4]],
		                    [rg_ten[3], rg_ten[1],rg_ten[5]],
		                    [rg_ten[4], rg_ten[5],rg_ten[2]]])

		#calculate parameters from the tensor of gyration
		
		L = np.linalg.eigvals(matrix)
		L.sort()
		L=L[::-1]
		Rg = (L[0]+L[1]+L[2])**0.5
		K1 = (L[1]+L[2])/(L[1]+L[0])
		K2 = (L[0]+L[2])/(L[0]+L[1])

		return [Rg, K1, K2, L[0], L[1], L[2]]

	def CN(
		self, 
		pos_A : np.ndarray | None = None, 
		pos_B : np.ndarray | None = None,
		dist : np.ndarray | None = None,
		d0 : float = 3.5,
		dmax  : float = 4.5
		) -> float | np.ndarray :
		"""
        Calculate the coordination number (CN) between atom group A and B,
        or within atom group A.
        
        Parameters
        ----------
        pos_A : np.ndarray of shape (N_A, 3) | None, optional
		    Either `pos_A` or `dist` must be available to perform calculations.
            Coordinates of the first atom group (Group A). 
            If provided, it overrides any `pos_A` saved in the class instance. 
            Either `pos_A` or `dist` must be available to performe calculatons.
        pos_B : np.ndarray of shape (N_B, 3) | None, optional
            Coordinates of the second atom group (Group B). 
            Required to calculate the CN between different groups. 
		    If provided, it overrides any `pos_B` saved in the class instance.
        dist : np.ndarray of shape (N_A, N_B) or (N_A, N_A) | None, optional
            Pre-calculated distance matrix. 
            Shape must be `(N_A, N_A)` for intra-group calculations (Group A only),
            or `(N_A, N_B)` for inter-group calculations. 
		    If provided, it overrides `pos_A` and `pos_B`.
        d0 : float, optional
            The inner cutoff parameter for the smoothing function.
        dmax : float, optional
            The outer cutoff parameter for the smoothing function.

        Raises
        ------
        ValueError
            Can't calculate CN without either pos_A/pos_B or dist provided.
        """
		pos_A, pos_B = self._set_pos(pos_A, pos_B)

		if all(v is None for v in [pos_A, pos_B, dist]):
			raise ValueError("All three pos_A, pos_B, and dist cannot be None.")

		if dist is None:
			dist = self.distance_matrix(pos_A, pos_B)
			
		else:
			dist = np.atleast_2d(dist)
			
		N_A = dist.shape[0]
		N_B = dist.shape[1]

		CN = np.zeros(N_A)

		for i in range(N_A):
			for j in range(N_B):

				if dist[i, j] == 0:
				    continue

				elif dist[i, j] <= d0:
				    CN[i] += 1

				elif d0 < dist[i, j] <= dmax:
				    CN[i] += self._cos(dist[i, j], dmax, d0)

		if len(CN)==1:
			return CN[0]

		else:
			return CN

	def rdf(
		self,
		pos_A : np.ndarray | None = None, 
		pos_B : np.ndarray | None = None,
		dist : np.ndarray | None = None,
		bins : int = 250,
		r_min : float = 0.0,
		r_max : float = 8.0
		) -> tuple[np.ndarray, np.ndarray]:
		"""
		Calculate the Radial Distribution Function (RDF) against distance (bin centers).
		
		Parameters
		----------
        pos_A : np.ndarray of shape (N_A, 3) | None, optional
            Coordinates of the first atom group (Group A). 
            If provided, it overrides any `pos_A` saved in the class instance. 
            Either `pos_A` or `dist` must be available to performe calculatons.
        pos_B : np.ndarray of shape (N_B, 3) | None, optional
            Coordinates of the second atom group (Group B). 
            Required to calculate the RDF between different groups. 
		    If provided, it overrides any `pos_B` saved in the class instance.
        dist : np.ndarray of shape (N_A, N_B) or (N_A, N_A) | None, optional
            Pre-calculated distance matrix. 
            Shape must be `(N_A, N_A)` for intra-group calculations (Group A only),
            or `(N_A, N_B)` for inter-group calculations. 
		    If provided, it overrides `pos_A` and `pos_B`.
		bins : int, optional
		    Number of bins used to discretize the distance range.
		r_min : float, optional
		    The lower limit of distance for the RDF calculation.
		r_max : float, optional
		    The upper limit of distance for the RDF calculation.
		
		Returns
		-------
		tuple[np.ndarray, np.ndarray]
		    - rdf : np.ndarray
		        The calculated RDF values.
		    - bin_centers : np.ndarray
		        The corresponding bin centers (in distance units).
		
		Raises
		------
		ValueError
		    Can't calculate RDF without either pos_A/pos_B or dist provided.
		"""
		pos_A, pos_B = self._set_pos(pos_A, pos_B)

		if all(v is None for v in [pos_A, pos_B, dist]):
			raise ValueError("All three pos_A, pos_B, and dist cannot be None.")

		if dist is None:
			dist = self.distance_matrix(pos_A, pos_B)
			dist[dist == 0] = np.nan
		
		N_A = dist.shape[0]

		bin_edges = np.linspace(r_min, r_max, bins + 1)
		bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
		dV = 4/3*np.pi*(bin_edges[1:]**3 - bin_edges[:-1]**3)

		hist = np.histogram(dist, bins=bin_edges)[0]
		rdf = hist/N_A/dV

		return rdf, bin_centers


	def x_min_distance(
		self,
		pos_A : np.ndarray | None = None, 
		pos_B : np.ndarray | None = None,
		dist : np.ndarray | None = None,
		x : int = 6 
		) -> float:
		"""
		Calculate xth minimal distance. By default, x=6.
		
		Parameters
		----------
		pos_A : np.ndarray of shape (N_A, 3) | None, optional
            Coordinates of the first atom group (Group A). 
            If provided, it overrides any `pos_A` saved in the class instance. 
            Either `pos_A` or `dist` must be available to performe calculatons.
        pos_B : np.ndarray of shape (N_B, 3) | None, optional
            Coordinates of the second atom group (Group B). 
            Required to calculate the xth minimal distance between different groups. 
		    If provided, it overrides any `pos_B` saved in the class instance.
        dist : np.ndarray of shape (N_A, N_B) or (N_A, N_A) | None, optional
            Pre-calculated distance matrix. 
            Shape must be `(N_A, N_A)` for intra-group calculations (Group A only),
            or `(N_A, N_B)` for inter-group calculations. 
		    If provided, it overrides `pos_A` and `pos_B`.
		x : int, optional
		    The number nearest neighbour
		    (1 -- the nearest, 2 -- second nearest etc.)
		
		Returns
		-------
		float
		    The parameter value in distance units.
		
		Raises
		------
		ValueError
		    Either pos_A/pos_B or dist must be provided.
		"""
		pos_A, pos_B = self._set_pos(pos_A, pos_B)

		if all(v is None for v in [pos_A, pos_B, dist]):
			raise ValueError("Either pos_A/pos_B or dist must be provided.")

		if dist is None:
			dist = self.distance_matrix(pos_A, pos_B)
			dist[dist == 0] = np.nan

		dist_sort = np.sort(dist, axis=1)
		x_min_distance = np.min(dist_sort[:, x-1])
		x_min_distance_at_id = np.where(dist_sort[:, 5]==x_min_distance)[0][0]

		return x_min_distance, x_min_distance_at_id


	def distance_matrix(
		self,
		pos_A : np.ndarray | None = None,
		pos_B : np.ndarray | None = None
		) -> np.ndarray:
		"""
		Calculate distance martix between atom group A and B (returns shape (N_A, N_B)),
		or within atom group A (returns shape (N_A, N_A))
		
		Parameters
		----------
		pos_A : np.ndarray of shape (N_A, 3) | None, optional
			Coordinates of the first atom group (Group A).
		    Must be provided if is not saved in the class instance.
            If provided, it overrides any `pos_A` saved in the class instance. 
        pos_B : np.ndarray of shape (N_B, 3) | None, optional
            Coordinates of the second atom group (Group B). 
            Required to calculate the distance matrix between different groups. 
		    If provided, it overrides any `pos_B` saved in the class instance.
		
		Returns
		-------
		np.ndarray
		    distance matrix.
		"""
		
		pos_A, pos_B = self._set_pos(pos_A, pos_B)

		if all(v is None for v in [pos_A, pos_B]):
			raise ValueError("At least one positional argument must be provided.")

		if pos_B is None:
			self_dist = distances.self_distance_array(pos_A, box=self.box)
			dist = self.self_dist_to_dist(self_dist, len(pos_A))

		else:
			dist = distances.distance_array(pos_A, pos_B, box=self.box)

		dist = np.atleast_2d(dist)


		return dist


	# ==========================================
	# HELPER METHODS (Internal)
	# ==========================================


	def _pbc(self, dist):
		"""
		Apply Periodic Boundary Conditions (PBC). Only orthogonal box:
		alpha, beta, gamma = 90, 90, 90. 
		
		Parameters
		----------
		dist : TYPE
		    non-pbc distance
		
		Returns
		-------
		TYPE
		    pbc distance
		"""
		dist = dist - self.box[:3] * np.rint(dist / self.box[:3])

		return dist

	@staticmethod
	def _cos(dist, dmax, d0):
	    """
	    Smoothing function for CN calculation.
	    
	    """
	    return 1/2*(np.cos((dist-d0)/(dmax-d0)*np.pi) + 1)

	@staticmethod
	def self_dist_to_dist(self_dist, N_A):
	    k = 0

	    dist = np.zeros((N_A, N_A))
	    for i in range(N_A):
	        for j in range(i + 1, N_A):
	            dist[i, j] = self_dist[k]
	            dist[j, i] = self_dist[k]
	            k += 1


	    return dist


	def _set_pos(self, pos_A, pos_B):

		if pos_A is None:
			pos_A = self.pos_A

		if pos_B is None:
			pos_B = self.pos_B


		return pos_A, pos_B
