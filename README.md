# pdbcop
Get details about a protein X-ray structure model from a PDB file.

- Print list of non standard residues in PDB file
- Print minimum, maximum, and average B factor values in chains and for the whole structure
- List all atoms with occupancy close to 0
- List all atoms with occupancy different from 0.00 and 1.00
- List all residues with alternate conformations
- Output sequence of protein chains in fasta format
- Compare two PDB files - usually before and after refinement - and print out otoms that differ in position and B factor more than a certain threshold

```
pdbcop [-l] [-d] [-s] file1 (file2) [shift_tol] [delta_B_tol]
```

-l long list: lists all atoms with occupancy different from 0.00 and 1.00  
-s prints out sequence  
-d compares atom positions and B-factors in two PDB files  
shift_tol: minimal shift to be printed, default value: 0.2  
delta_B_tol: minimal difference in B to be printed, default value: 3.5 

*Disclaimer: No one has touched the code since 2015.*
