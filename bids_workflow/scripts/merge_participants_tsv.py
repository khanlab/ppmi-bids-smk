import pandas as pd


df_list = [pd.read_table(tsv) for tsv in snakemake.input.tsv_files]

df = pd.concat(df_list)
print(df)
df.to_csv(snakemake.output.tsv_file,sep='\t',index=False)

