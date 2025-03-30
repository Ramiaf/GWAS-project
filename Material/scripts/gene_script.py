import pysam
import argparse

def find_genes_with_snp(tabix_file, chrom, pos):
    if chrom == 23:
        chrom = 'X'
    region = f"{chrom}:{pos-1}-{pos}"     #pos-1 and pos to narrow down the search as much as possible
    genes_containing_snp = []
    
    for entry in tabix_file.fetch(region=region):     #it will retrieve all entries/genes that overlap with the specified region
        fields = entry.split('\t')
        gene_chrom = fields[0]
        gene_start = int(fields[1])
        gene_end = int(fields[2])
        gene_name = fields[3]
        gene_NCBI_id = fields[5]
        genes_containing_snp.append((gene_name, gene_NCBI_id, gene_start, gene_end))
    if len(genes_containing_snp) == 0:
        return "None"
    if len(genes_containing_snp) == 1:
        return genes_containing_snp[0]
    else:    #case of nested genes
        gene_of_interest = genes_containing_snp[0]
        min_start = gene_of_interest[2]
        min_end = gene_of_interest[3]
        for gene in genes_containing_snp:
            if gene[2] < min_start and gene[3] < min_end:
                gene_of_interest = gene
                min_start = gene[2]
                min_end = gene[3]
        return gene_of_interest
    
compressed_file="HS_genes/HS_final_sorted_genes_copy.txt.gz"
tabix_file = pysam.TabixFile(compressed_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Find genes containing SNP.')
    parser.add_argument('--snp_chrom', type=int, required=True, help='Chromosome of the SNP')
    parser.add_argument('--snp_pos', type=int, required=True, help='Position of the SNP')

    args = parser.parse_args()

    genes = find_genes_with_snp(tabix_file, args.snp_chrom, args.snp_pos)
    print(genes)
