import numpy as np
from infretis.classes.orderparameter import OrderParameter
import MDAnalysis as mda
from ionicluster import StructureAnalyser
from ionicluster import ClusterAnalyser

u = mda.Universe('lammps_input/lammps.data', format='DATA')

Ca = u.select_atoms('type 1 ')
ca_id_list = Ca.indices 

C = u.select_atoms('type 2')
c_id_list = C.indices 

Ow = u.select_atoms('type 4')
ow_id_list = Ow.indices

struct = StructureAnalyser(init = True)
cl = ClusterAnalyser(init = True)


class OrderX(OrderParameter):

    def __init__(self):
        super().__init__()

    def calculate(self, system):

        # set up system variables
        pos = system.pos.astype(np.float32)
        box = np.array([*system.box,90,90,90]).astype(np.float32)

        ca_pos = pos[ca_id_list]
        c_pos = pos[c_id_list]
        ow_pos = pos[ow_id_list]

        struct.update(box=box)

        ca_c_matr = struct.distance_matrix(ca_pos, c_pos)
        op, central_at_id = struct.x_min_distance(dist=ca_c_matr) # Ca-C 6th min OP

        # calculate CN

        cn_cac_c = struct.CN(dist=ca_c_matr[central_at_id], d0=3.5, dmax=4.5)  # Ca cental -- C CN
        central_at_pos = ca_pos[central_at_id]
        cn_cac_ca = struct.CN(pos_A=central_at_pos, pos_B=ca_pos, d0=4.5, dmax=5.5) # Ca cental -- Ca CN
        cn_cac_ow = struct.CN(pos_A=central_at_pos, pos_B=ow_pos, d0=3.5, dmax=5.5) # Ca cental -- Ow CN

        # cluster analysis

        cl.update(box, ca_pos, c_pos)
        n_cl,_ = cl.count()
        CI_cl_at_pos, _, CI_cl_n_at = cl.CI_cluster(central_at_id=central_at_id) # CI_cl_n_at -- cluster size

         # shape analysis

        Rg, K1, K2, _, _, _ = struct.Rg_K1_K2(pos=CI_cl_at_pos)
        
        return [-1*op, central_at_id,
                cn_cac_c, cn_cac_ca, cn_cac_c+cn_cac_ca, cn_cac_ow,
                CI_cl_n_at, n_cl, Rg, K1, K2]