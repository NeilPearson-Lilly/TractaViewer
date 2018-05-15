# This is a script to make BLAST output into the identity score and likelihood matrices we want to try using.
# This matrix will be fairly hefty, but still preferable to actually doing alignments on-the-fly.

import pandas as pd
import re
import json

res = pd.read_table("results.tsv", names=["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart",
                                          "qend", "sstart", "send", "evalue", "bitscore"])

ident_matrix = {}
eval_matrix = {}
c = 0
for i, row in res.iterrows():
    # So we know how quick things are going...
    c += 1
    if c % 1000000 == 0:
        print(c)
    qry_id = re.split("[\|-]", row['qseqid'])[1]
    hit_id = re.split("[\|-]", row['sseqid'])[1]
    if qry_id not in ident_matrix:
        ident_matrix[qry_id] = {}
    if hit_id in ident_matrix[qry_id]:
        if row['pident'] > ident_matrix[qry_id][hit_id]:
            ident_matrix[qry_id][hit_id] = row['pident']
    else:
        ident_matrix[qry_id][hit_id] = row['pident']
    
    if qry_id not in eval_matrix:
        eval_matrix[qry_id] = {}
    if hit_id in eval_matrix[qry_id]:
        if row['evalue'] < eval_matrix[qry_id][hit_id]:
            eval_matrix[qry_id][hit_id] = row['evalue']
    else:
        eval_matrix[qry_id][hit_id] = row['evalue']

# idents = pd.DataFrame(ident_matrix)
# evals = pd.DataFrame(eval_matrix)
# 
# idents.to_csv("ident_matrix.csv")
# evals.to_csv("eval_matrix.csv")

# Still too big. Let's try JSON. 
with open("ident_matrix.json", 'w') as f:
    json.dump(ident_matrix, f)
with open("eval_matrix.json", 'w') as f:
    json.dump(eval_matrix, f)

