Overview of tree sequence recording
--------------------------------------------------------------

fwdpy11 0.2.0 added support for tree sequence recording (TSR), which is a method to both speed up
forward-time simulations and to track the entire genealogical history of the simulation.  The key reference
describing TSR is::

    Kelleher, Jerome, Kevin Thornton, Jaime Ashander, and Peter Ralph. 2018.
    Efficient Pedigree Recording for Fast Population Genetics Simulation.
    bioRxiv. https://doi.org/10.1101/248500.

The methods described in the 2018 paper rely heavily on concepts described in the preceding work::

    Kelleher, Jerome, Alison M. Etheridge, and Gilean McVean. 2016.
    Efficient Coalescent Simulation and Genealogical Analysis for Large Sample Sizes.
    PLoS Computational Biology 12 (5): e1004842.
