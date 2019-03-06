import pandas as pd

# This is a testing ground to read this staggeringly poorly formatted file.
lines = []
with open("Data/pdbtosp.txt") as file:
    for line in file: 
        # line = line.strip()
        lines.append(line)

prevcode = None
prevmethod = None
prevres = None
good_lines = []
for line in lines[24:-5]:
    pdb_code = line[0:5].strip()
    if not pdb_code:
        pdb_code = prevcode
    str_method = line[6:15].strip()
    if not str_method:
        str_method = prevmethod
    res = line[16:27].strip()
    if not res:
        res = prevres
    sp1 = line[28:39].strip()
    up1 = line[41:47].strip()
    l1 = [pdb_code, str_method, res, sp1, up1]
    # print(l1)
    good_lines.append(l1)
    if len(line) > 53:
        sp2 = line[53:65].strip()
        up2 = line[67:73].strip()
        l2 = l1 = [pdb_code, str_method, res, sp2, up2]
        # print(l2)
        good_lines.append(l2)
    
    prevcode = pdb_code
    prevmethod = str_method
    prevres = res
    
df = pd.DataFrame(good_lines, columns=['PDB ID', 'Method', 'Resolution', 'SwissProt ID', 'UniProt ID'])
print(df)
df.to_csv("Data/pdb_to_uniprot.csv", index=False)
