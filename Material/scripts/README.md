Both `gene_script.py` and `vcf_script.py` are used to retrieve SNP information for the raw SNPs of the dataframe.<br><br>
However, `gene_script.py` is used after `vcf_script.py` is complete because `vcf_script.py` only returns SNP info for known SNPs in the NCBI database.<br><br>
`gene_script.py` will search for the SNPs that were not identified by `vcf_script.py`, find their corresponding genes, and then extract SNP information.<br><br>
(this allows for discovery of novel SNPs that could be highly associated with the trait of interest).
