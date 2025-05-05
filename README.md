### AOA_urea

This repository hosts the scripts associated with the sequence analyses for the work described in Stuehrenberg, Kitzinger, et al. 

#### scripts
The directory "assembly refine" contains a bash script for iterative read mapping on a metagenomic bin and subsequent assembly of the mapping reads with the goal to increase bin completeness and contiguity.

The directory "read searching" contains a python script to perform marker gene based analyis on unassembled read datasets. This script can be used with the databases in the directory "marker gene databases".

The directory "marker gene databases" contains a bash script, and three helper scripts in perl, to build marker gene protein databases. The subdirectory "databases" contains databases for the amoA/pmoA, ureC, and dur3 genes built using these scripts and the total protein complement of a dereplicated genome set of GTDB release 207 (https://gtdb.ecogenomic.org/) and the GEM OTU set (https://www.nature.com/articles/s41587-020-0718-6).
