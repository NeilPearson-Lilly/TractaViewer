import pandas as pd
import biomart
from io import StringIO

# A quick script to make this dataset. By default, it identifies genes by their ChEMBL ID, which is not much good to
# me. Replace it with UniProt ID or something - anything - else. 

orig = pd.read_table("Data/chembl_drugtargets-18_13_46_02.txt")

server = biomart.BiomartServer("http://uswest.ensembl.org/biomart")
dataset = server.datasets['hsapiens_gene_ensembl']
# dataset.show_filters()
# dataset.show_attributes()

chembl_ids = orig['TARGET_CHEMBL_ID'].dropna().unique().tolist()



attrs = ['chembl', 'hgnc_symbol',
         # 'ensembl_gene_id', 'uniprot_gn'
         ]
response = dataset.search({'filters': {'chembl': chembl_ids[:500]}, 
                           'attributes': attrs})
data = "\n".join([i.decode('utf-8') for i in response.iter_lines()])
df = pd.read_table(StringIO(data), names=attrs)

response = dataset.search({'filters': {'chembl': chembl_ids[500:]}, 
                           'attributes': attrs})
data = "\n".join([i.decode('utf-8') for i in response.iter_lines()])
df2 = pd.read_table(StringIO(data), names=attrs)

df3 = pd.concat([df, df2], ignore_index=True)
df3.columns = ["TARGET_CHEMBL_ID", "HGNC Name"]
print(df3)

final = orig.merge(df3, how='left', left_on="TARGET_CHEMBL_ID", right_on="chembl")

final.to_csv("Data/chembl_drugtargets_named-18_13_46_02.csv")

