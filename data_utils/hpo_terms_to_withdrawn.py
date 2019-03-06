import pandas as pd

# This produces a lookup table fetching the number of genes and number of withdrawn drugs falling into these categories.
# The values from WITHDRAWN aren't easy to get at from a script, so I'm hardcoding them. WITHDRAWN doesn't update very 
# often anyway.

hpo_terms = ['Hepatic', 'Cardiovascular', 'Hematological', 'Neurological', 'Dermatological', 'Cancer', 'Renal', 
             'Respiratory', 'Reproductive', 'Opthalmic', 'Gastrointestinal', 'Muscular']
hpo_ids = ['HP:0001392', 'HP:0001626', 'HP:0001871', 'HP:0000707', 'HP:0000951', 'HP:0002664', 'HP:0000077', 
           'HP:0002086', 'HP:0000119', 'HP:0000478', 'HP:0025031', 'HP:0003011']
withdrawn_classes = ['Abnormality of the liver', 'Abnormality of the cardiovascular system', 
                     'Abnormality of blood and blood-forming tissues', 'Abnormality of the nervous system', 
                     'Abnormality of the skin', 'Neoplasm', 'Abnormality of the kidney', 
                     'Abnormality of the respiratory system', 'Abnormality of the genitourinary system', 
                     'Abnormality of the eye', 'Abnormality of the digestive system', 'Abnormality of the musculature']
# These numbers come from this figure: http://cheminfo.charite.de/withdrawn/d3/tox_type_vs_atc.html
# There is no way to automatically extract these, as far as I can tell.  
drug_tox_counts = [60, 48, 31, 27, 24, 22, 12, 9, 9, 6, 5, 4] 


wd = pd.read_csv("Data/withdrawn_withdrawn_compounds.csv")
hp = pd.read_table("Data/HPO_genes_to_diseases.txt", skiprows=[0], names=['entrez', 'gene', 'disease'])

num_genes = []
for i, j in zip(hpo_terms, hpo_ids):
    associated_diseases = pd.read_csv("Data/HPO_assocs/" + i + "_diseases.csv", 
                                      skiprows=[0], 
                                      names=['disease_id', 'disease_name'])
    genes = list(set(hp[hp['disease'].isin(associated_diseases['disease_id'])]['gene'].tolist()))
    num_genes.append(len(genes))


df = pd.DataFrame({'Term name': hpo_terms,
                   'HPO annotation': withdrawn_classes,
                   'HPO ID': hpo_ids,
                   'Num. genes': num_genes,
                   'Withdrawn drug count': drug_tox_counts})
print(df)

