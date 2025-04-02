# This folder contains genes files (unsorted and sorted) retrieved from Biomart using python
  - HS_genes.txt contains all the genes extracted from Biomart (unsorted).
  - HS_initial_sorted_genes.txt contains the sorted genes using pandas (not enough since non-numerical genes and haplotypes are still not sorted properly).
  - This file is then split into 2 files:
      - HS_sorted_CME.txt contains all the numerical chromosomes in a sorted fashion.
      - HS_unsorted_haplotypes.txt contains the remaining non-numerical genes (unsorted).
  - HS_sorted_haplotypes.txt contains the sorted haplotypes.
  - HS_final_sorted_genes.txt is the result of merging "HS_sorted_haplotypes.txt" and "HS_sorted_CME.txt".<br>
#### Now we can index this file and use it to quickly extract gene and SNP information.
