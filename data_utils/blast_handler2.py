# This is a script to make BLAST output into the identity score and likelihood matrices we want to try using.
# This matrix will be fairly hefty, but still preferable to actually doing alignments on-the-fly.

import pandas as pd
import re
import json

res = pd.read_table("results.tsv", names=["qseqid", "sseqid", "pident", "qlen", "slen", "length", "nident", "mismatch", 
                                          "gapopen", "evalue", "bitscore"])

ident_matrix = {}
# qlens = {}
eval_matrix = {}
c = 0
for i, row in res.iterrows():
    # So we know how quick things are going...
    c += 1
    if c % 1000000 == 0:
        print(c)
    qid = re.split("[\|-]", row['qseqid'])[1]
    sid = re.split("[\|-]", row['sseqid'])[1]
    
    # qlens[qid] = row["qlen"]
    if row["length"] >= (row["qlen"] / 10) * 8:
        ident = (row["nident"] / row["qlen"]) * 100
        ident_inverse = (row["nident"] / row["slen"]) * 100
        if qid not in ident_matrix:
            ident_matrix[qid] = {}
        if sid in ident_matrix[qid]:
            if ident > ident_matrix[qid][sid]:
                ident_matrix[qid][sid] = ident
        else:
            ident_matrix[qid][sid] = ident
        # Do that the other way round too
        if sid not in ident_matrix:
            ident_matrix[sid] = {}
        if qid in ident_matrix[sid]:
            if ident > ident_matrix[sid][qid]:
                ident_matrix[sid][qid] = ident_inverse
        else:
            ident_matrix[sid][qid] = ident_inverse
        
        # And do e-values
        if qid not in eval_matrix:
            eval_matrix[qid] = {}
        if sid in eval_matrix[qid]:
            if row['evalue'] < eval_matrix[qid][sid]:
                eval_matrix[qid][sid] = row['evalue']
        else:
            eval_matrix[qid][sid] = row['evalue']
    

# for qid in ident_matrix:
#     for sid in ident_matrix[qid]:
#         ident_matrix[qid][sid] = (ident_matrix[qid][sid] / qlens[qid]) * 100

# def get_uniprot_id_back(id):
#     return re.split("[\|-]", id)[1]
# 
# 
# # Remove some useless formatting
# res['qid'] = res["qseqid"].apply(get_uniprot_id_back)
# res['sid'] = res["sseqid"].apply(get_uniprot_id_back)
# 
# ident_matrix = {}
# eval_matrix = {}
# c = 0
# for qid in res["qid"].unique():
#     c += 1
#     if c % 100 == 0:
#         print(str(c) + " of " + str(len(res["qid"].unique())))
#     for sid in res["sid"].unique():
#         rel = res[(res["qid"] == qid) & (res["sid"] == sid)]
#         if len(rel) > 0:
#             ident = (rel['nident'].sum() / float(rel['qlen'].values[0])) * 100
#             try:
#                 ident_matrix[qid][sid] = ident
#             except KeyError:
#                 ident_matrix[qid] = dict()
#                 ident_matrix[qid][sid] = ident
#             # Not sure this is the best way to do this, but let's go with it for now.
#             e = rel['evalue'].min()
#             try:
#                 eval_matrix[qid][sid] = e
#             except KeyError:
#                 eval_matrix[qid] = dict()
#                 eval_matrix[qid][sid] = e

# Still too big. Let's try JSON. 
with open("ident_matrix.json", 'w') as f:
    json.dump(ident_matrix, f)
with open("eval_matrix.json", 'w') as f:
    json.dump(eval_matrix, f)

c1 = 0
c2 = 0
for i in ident_matrix:
    for j in ident_matrix[i]:
        c1 += 1
        if ident_matrix.get(j).get(i):
            c2 += 1
print("Total number of records in matrix: " + str(c1))
print("Number of records with a symmetrical figure: " + str(c2))
