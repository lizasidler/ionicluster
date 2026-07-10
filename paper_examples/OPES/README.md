# Paper Examples

This directory contains scripts designed for the post-processing of simulation (MD) trajectories (OPES in this case, but adaptable to any set of chronologically numbered trajectories) for calcium carbonate simulations in water.

## Scripts

*   **`nb_list_test.py`**  
    Demonstrates how different kinds of the Coordination Number (Ca--Ca CN, C--Ca CN, ...) can be calculated for all selected atoms (all Ca and all C in this case) across all trajectories using `StructureAnalyser`.

*   **`rdf_test.py`**  
    Demonstrates how to calculate different cluster Radial Distribution Functions (RDFs) and Pair Distribution Functions (PDFs) across all trajectories using `StructureAnalyser` and `ClusterAnalyser`.