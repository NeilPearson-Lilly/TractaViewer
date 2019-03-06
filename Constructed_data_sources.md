# Constructed data sources

A small number of data sources used by TractaViewer are constructed by scripts provided in the `data_utils/` directory. Most data sources bundled with TractaViewer can be updated, should the need arise, by simply replacing the existing files with newer versions - provided that the source in question maintains a consistent format between versions. However, a small number of sources require more extensive modification before they can be used.

### PDB -> UniProt lookup table

In order to quickly look up associated UniProt IDs from PDB IDs, we need to use the script `data_utils/pdb_to_uniprot.py` to create a lookup table. The input to this script is [this file](https://www.uniprot.org/docs/pdbtosp.txt) from UniProt. Output is `Data/pdb_to_uniprot.csv`. 

### Human Phenotype Ontology -> withdrawn drugs lookup table

We produced a table linking broad categories of human disease phenotypes to affected tissue categories in the WITHDRAWN database. This is the output of the script `data_utils/hpo_terms_to_wihdrawn.py`. Inputs are the top-level toxicity type classifications used by the WITHDRAWN database which can be viewed [here](http://cheminfo.charite.de/withdrawn/d3/tox_type_vs_atc.html), and the most closely corresponding high-level terms in the Human Phenotype Ontology. The relevant data from the HPO is bundled with TractaViewer, and is linked in the data manifest. 

Since this information is not conveniently accessible in an automated manner from WITHDRAWN, and since the final lookup table produced is not large, this table was produced primarily by hand.

Output is `Data/HPO_termname_lookup.csv`. 

### Genome-wide identity and e-value matrices

To compute similarity of protein sequences to other proteins across the whole protein-coding genome, we ran an all-vs-all BLASTp using all protein sequences from UniProt, with the following parameters:

- outfmt "6 qseqid sseqid pident qlen slen length nident mismatch gapopen evalue bitscore"

Output from this alignment is read by the script `data_utils/blast_handler2.py`, which produces the files `Data/ident_matrix.json` and `Data/eval_matrix.json`.

### Withdrawn drug information

A list of withdrawn drugs is available from the [WITHDRAWN database](http://cheminfo.charite.de/withdrawn/links.html). It comes in the form of an SDF file, which can be converted to a CSV file using `data_utils/sdf_to_csv.py`

### Named ChEMBL drug targets

ChEMBL maintains a list of drugs and their intended target genes, available from its [download page](https://www.ebi.ac.uk/chembl/drug/targets). In this file, genes are identified only by their ChEMBL IDs. We created `data_utils/make_chembl_targets_data.py` to fill in other IDs from Ensembl's BioMart. 
