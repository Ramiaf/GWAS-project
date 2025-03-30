#Script to create a GWAS dataframe with all data I retrieved for the raw SNPs
import pysam
import pandas as pd

def find_snp(vcf_file):                                    
    for i in range(len(GWAS_df)):                              
        chrom = GWAS_df['CHROM'][i]                            
        SNP_pos = GWAS_df['GENPOS'][i]
        records = vcf_file.fetch(chrom, SNP_pos-1, SNP_pos)
        #example of record:
        #CHROM     POS        ID    REF   ALT   QUAL  FILTER  INFO
        #1   55550  rs1234567      G     A     .      .     'dict object'

        for record in records:   #loop needed because 'records' is an iterator object (even if only one match)
            record_info = dict(record.info)  #example entry{'RS': 775809821, 'RSPOS': 10020, 'dbSNPBuildID': 144, 'SSR': 0, 'SAO': 0, 'VP': '0x050000020005000002000200', 'GENEINFO': 'DDX11L1:100287102', 'WGT': 1, 'VC': 'DIV', 'R5': True, 'ASP': True}
            if 'GENEINFO' in record_info:
                gene_info = record_info['GENEINFO'].split(":")
                GWAS_df.loc[(GWAS_df["CHROM"] == record.chrom) & (GWAS_df['GENPOS'] == record.pos), "Gene_name"] = gene_info[0]
                GWAS_df.loc[(GWAS_df["CHROM"] == record.chrom) & (GWAS_df['GENPOS'] == record.pos), "NCBI_ID"] = gene_info[1]
            GWAS_df.loc[(GWAS_df["CHROM"] == record.chrom) & (GWAS_df['GENPOS'] == record.pos), "ID"] = record.id

def search_gene(CHROM, SNP):
    # Whether the Gene name is known or not, we still need to go back to the 'genes file' to retrieve the start and end positions of this gene
    command = f'wsl python3 scripts/gene_script.py --snp_chrom {CHROM} --snp_pos {SNP}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    # Check the output
    out = result.stdout.strip()
    if out == "None":
        return "None"
    else:
        try:
            result_list = ast.literal_eval(out)  # Safely evaluate the string to a list
        except (SyntaxError, ValueError) as e:
            print(f"Error parsing result_str: {result_str}")
        return result_list


GWAS_df = pd.read_csv("regenie_output/GWAS_results.txt", sep='\t')

GWAS_df["Gene_name"] = 'NA'  #Use this column for when you find the Gene of the SNP; replace NA with gene name.
GWAS_df['NCBI_ID'] = 'NA'
GWAS_df["GENPOS"] = GWAS_df["GENPOS"].astype('int64')

GWAS_df["CHROM"] = GWAS_df["CHROM"].replace({"X":"23"})
GWAS_df["CHROM"] = GWAS_df["CHROM"].astype('int64')

GWAS_df["color"] = GWAS_df['CHROM'].apply(lambda x : 'blue' if x%2==0 else 'black')   #For coloring purposes of manhattan plot
GWAS_df['def_color'] = GWAS_df['color']


compressed_file="00-All.vcf.gz"
vcf_file = pysam.VariantFile(compressed_file)
find_snp(vcf_file)

#I am going to loop over the entire GWAS_df ONCE, find the gene_names of the unidentified SNPs using the search_gene() function.
for index, row in GWAS_df.iterrows():
    if type(row['Gene_name']) == float:  #if it is null object (no gene name)
        gene_info = search_gene( row['CHROM'] , row['GENPOS'] )
        
        if gene_info != 'None':
            GWAS_df.loc[index, 'Gene_name'] = gene_info[0]
            if gene_info[1] != '':
                GWAS_df.loc[index, 'NCBI_ID'] = int(float(gene_info[1]))


#We want to add a column that will keep track of the cumulative positions of the SNPs for manhattan plot
start_pos_relative_to_genome = 0
cumulative_list = []

#groupby("CHROM") groups the CHROM column where the name of the group is the chromosome and the group is a subset of the dataframe whose chrom is of that group
for chrom, group in GWAS_df.groupby("CHROM"):
    cumulative_list.append(group["GENPOS"] + start_pos_relative_to_genome)    #store the positions of the SNPs relative to the genome, starts with 0 then the max SNP of chrom 1 and so on..
    start_pos_relative_to_genome += group["GENPOS"].max()
GWAS_df["cumulative_pos"] = pd.concat(cumulative_list)

GWAS_df["CHROM"] = GWAS_df["CHROM"].astype(str)
GWAS_df["CHROM"] = GWAS_df["CHROM"].replace({"23":"X"})

GWAS_df.to_csv('GWAS_df_updated.csv', index=False)
