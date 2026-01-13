# Mountain Pine Beetle Spatial Prioritisation of Slow-the-Spread

![Results plot](Plots/1-Base/Fig4.eps)

This work accompanies Hudgins et al. (2026) *The value of looking ahead: comparing conventional and strategic Mountain Pine Beetle (Dendroctonus ponderosae) management policies in North America*

This code was written by Emma Hudgins and Denys Yemshanov

R scripts are numbered and should be run in order.

Results folders and Plots folders are numbered by the indexed sensitivity analysis scenarios (where 1 is the best guess of parameter values).

The GAMS script (launched using the bash file between R scripts 1 and 2) requires a GAMS software license (https://www.gams.com) and is shared for transparency of model formulation in spite of the non-open source nature of this software. 
Ths script can be converted into an open source MILP formulation and run freely, but note that solving times will be much higher.

Some raw data files are not included because they contain private government surveillance information.

Files are described in comments within the scripts included.
