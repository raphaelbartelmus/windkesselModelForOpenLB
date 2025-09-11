# windkesselModelForOpenLB
This repository contains an implementation of the three element windkessel (3WK) model that can be used in OpenLB.

# Usage:
1. Download [OpenLB 1.7](https://doi.org/10.5281/zenodo.10684609)
2. Copy *`wkfct3d.h`* into *`/olb-1.7r0/src/`*
3. Create a new directory in one of the example folders e.g. *`/olb-1.7r0/examples/turbulence/windkessel3d/`* and copy *`windkessel3d.cpp`* and *`Makefile`* into this directory
4. The code expects an *`.stl`* file with the geometry to be in the same directory. If you possess the original geometry, simply copy it into the folder. If not, you must:
    - change the name of the .stl file in line 465 to match the name of your file
    - adjust the definition of the location inlets and outlets **twice**, in **`prepareGeometry(...)`** and ll. 507-562
    - adjust the definition of inlets and outlets in **`prepareLattice(...)`**
    - adjust the number of 3WK outlets in the definition and use of **`rkUpdtater(...)`** and in **`setBoundatryValues(...)`**
    - adjust the 3WK-parameters in ll. 87-114, 447-456 and in the definition of **`rkUpdtater(...)`**
5. Set the parameters in ll. 33-67 to your needs
6. compile the program with `make` and run it
