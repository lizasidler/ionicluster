# Paper Examples

This directory contains scripts designed for the post-processing of infRETIS paths for calcium carbonate simulations in water.

## 🧪 Test System & Example Paths

The scripts are pre-configured to be tested with a system defined by $\chi_{\rm ion/w} = 0.01$ and $N_{\rm ion} = 72$, the Ca–C order parameter (OP). 

The exemplary path data is located in the `chi_0.01_nion_72_cac_OP/` folder:
*   **`chi_0.01_nion_72_cac_OP/load/0`**: Trajectory data for the [0-] path.
*   **`chi_0.01_nion_72_cac_OP/load/1`**: Trajectory data for the reactive path.

### ⚠️ Optimizing File Sizes (`cut_water.py`)
The raw trajectory files are relatively large (~400 MB). In order to optimize storage and significantly speed up your local test runs, water moleculed were removed from the trajectories given in `chi_0.01_nion_72_cac_OP/` using cut_water.py.

## Scripts

*   **`nb_list_test.py`**  
    Demonstrates how different kinds of the Coordination Number (Ca--Ca CN, C--Ca CN, ...) can be calculated for all selected atoms (all Ca and all C in this case) across all trajectories using `StructureAnalyser`.

*   **`rdf_test.py`**  
    Demonstrates how to calculate different cluster Radial Distribution Functions (RDFs) and Pair Distribution Functions (PDFs) across all trajectories using `StructureAnalyser` and `ClusterAnalyser`.

*   **`rg_hist_cl_size_test.py`**  
    Demonstrates how to calculate the Radius of Gyration ($R_g$) distribution per cluster size across all trajectories using `StructureAnalyser` and `ClusterAnalyser`.

*   **`cut_water.py`** 
    Demonstrates how to remove water molecules from the result trajectories in `load`. Allows for the great reduction in disk space.