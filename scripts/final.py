import pandas as pd

pheno_df = pd.read_excel("selected_data.xlsx", header=0, index_col=0)
geno_df = pd.read_csv("final/bacteria.prune.plink.raw", sep='\s+', header=0, index_col=0)
geno_df = geno_df.drop(columns=['IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'])
pheno_df.join(geno_df, how='inner').to_csv("final_data.csv")
