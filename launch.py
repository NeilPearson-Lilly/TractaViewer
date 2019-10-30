from glob import glob
from json import load, loads, dump
from os import remove, getcwd, mkdir
from os.path import expanduser, basename, exists, isfile, isdir
from re import sub, search, match, split, IGNORECASE
from shutil import move
from sys import argv, exit
from traceback import print_exc
from xml.parsers.expat import ExpatError
from collections import OrderedDict
from io import StringIO
from pprint import pprint
from time import sleep
from functools import reduce

import numpy.core._multiarray_umath
from numpy import nan
from pandas import set_option, read_excel, read_csv, read_json, read_html, read_table, to_numeric, DataFrame, \
    isnull, concat
from pweave import weave
from requests import get, Session
import xmltodict
import PyQt5
from PyQt5 import QtCore, QtWidgets
from PyQt5.QtCore import QThread, pyqtSignal, Qt, QItemSelectionModel, QVariant
from PyQt5.QtGui import QPixmap, QIcon, QIntValidator
from PyQt5.QtWidgets import QApplication, QWidget, QMainWindow, QDialog, QFileDialog, QMessageBox, QTableView, \
    QSplashScreen
from biomart import BiomartServer
from bs4 import BeautifulSoup
# import PyEntrezId.Conversion
from intermine.webservice import Service
from opentargets import OpenTargetsClient
from pandas import ExcelWriter
from pypdb import describe_pdb
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry
from sklearn.preprocessing import StandardScaler

import about
import bucketing_help
import disease_profile_dialog
import gui
import makereports_dialog
import settings
import textinsert

# set_option('display.max_columns', 100)
# set_option('display.width', 220)
set_option('display.max_columns', None)
set_option('display.max_rows', None)

# Useful variables
disease_association_score_threshold = 0.5
sheet_order = ['External data input', 'Basic information', 'Buckets', 'Bucket guide', 'GTEX', 'HPA expression', 
               'Barres mouse', 'Barres human', 'HPA Enrichment', 'Risk factors', 'Pathways', 
               'Rare disease associations', 'Common disease associations', 'MGI mouse', 'CanSAR', 'Protein structures', 
               'Antitargets', 'Core fitness', 'SM Druggability', 'AB-ability', 'Feasibility', 'Existing drugs', 
               'Drug toxicology', 'Pharos', 'GWAS', 'Protein-protein interactions', 'Literature']

# We're going to put these in a separate file and not display these in the GUI, because they bring everything to a crawl
# when we're working with more than a few targets. Worth pulling down though. 
unacceptably_large_sheets = ['Protein-protein interactions', 'Literature']

# The following columns should be represented as booleans in model views:
bool_cols = ['GPCR', 'Predicted membrane protein', 'Predicted secreted protein', 'Surfaceome membership (human)',
             'Surfaceome membership (mouse)', 'Mutational cancer driver genes',
             'COSMIC somatic mutations in cancer genes', 'Core fitness gene', 'CRISPR-screened core fitness gene',
             "Membrane proteins predicted by MDM", "MEMSAT3 predicted membrane proteins",
             "MEMSAT-SVM predicted membrane proteins", "Phobius predicted membrane proteins",
             "SCAMPI predicted membrane proteins", "SPOCTOPUS predicted membrane proteins",
             "THUMBUP predicted membrane proteins", "TMHMM predicted membrane proteins",
             "GPCRHMM predicted membrane proteins", "Secreted proteins predicted by MDSEC",
             "SignalP predicted secreted proteins", "Phobius predicted secreted proteins",
             "SPOCTOPUS predicted secreted proteins"]

# One from MGI: these are the phenotype classes one up from the root. 
high_level_phenotypes = ['adipose tissue', 'behavior/neurological', 'cardiovascular system', 'cellular',
                         'craniofacial', 'digestive/alimentary', 'embryo', 'endocrine/exocrine glands',
                         'growth/size/body', 'hearing/vestibular/ear', 'hematopoietic system', 'homeostasis/metabolism',
                         'integument', 'immune system', 'limbs/digits/tail', 'liver/biliary system', 'mortality/aging',
                         'muscle', 'nervous system', 'pigmentation', 'renal/urinary system', 'reproductive system',
                         'respiratory system', 'skeleton', 'taste/olfaction', 'neoplasm', 'vision/eye']


# Useful funcs
def uniq(lst):
    return list(set(lst))


def flatten(ls):
    # Recursively turn a list of lists into a simple list
    if isinstance(ls, list):
        if len(ls) == 0:
            return []
        first, rest = ls[0], ls[1:]
        return flatten(first) + flatten(rest)
    else:
        return [ls]


# We need a couple of lists available from the start. These can be accessed by looking through some folders for files.
hpa_toxfiles = glob("Data/HPA_tox_lists/*.csv")
hpa_all_tissues = []
for fl in hpa_toxfiles:
    tissue, datatype = [sub('\+', ' ', i).capitalize() for i in split('[_\.]', basename(fl))[0:2]]
    hpa_all_tissues.append(tissue)
hpa_all_tissues = sorted(uniq(hpa_all_tissues))

hpo_annotation_files = glob("Data/HPO_assocs/*_diseases.csv")
hpo_annotation_names = sorted(uniq([sub("_diseases.csv", '', basename(i)) for i in hpo_annotation_files]))


def deep_get(dictionary, keys, default=None):
    return reduce(lambda d, key: d.get(key, None) if isinstance(d, dict) else default, keys, dictionary)


class Downloader(QThread):
    progbar_update = pyqtSignal(int)
    status = pyqtSignal(str)
    got_data = pyqtSignal(object)
    warnings = pyqtSignal(str)
    
    def __init__(self, api_keys, literature, ppi, organism="human", gene_family_query_level="node"):
        super(QThread, self).__init__()
        self.api_keys = api_keys
        # These get tagged as True if the appropriate checkboxes are tagged in the GUI. We do this to give the option of
        # not downloading these large but potentially useful data-dumps. 
        self.get_literature = literature
        self.get_ppi = ppi
        
        self.disease_profile = None
        self.input_data = None
        self.organism = organism
        self.gene_family_query_level = gene_family_query_level
        # We have a number of hardcoded lists - fill them up in a sub so we can hide it for readability.
        self.toplevel = []
        self.secondlevel = []
        self.thirdlevel = []
        self.poss_existing_cols = []
        
        self.gtex_tissues = []
        self.barreslab_mouse_tissues = []
        self.barreslab_mouse_orig_cols = []
        self.barreslab_human_tissues = []
        self.general_tissues = []
        self.cell_types = []
        self.significance_threshold = 0.05
        self.diff_expression_threshold = 2
        # These come from http://juniper.health.unm.edu/tcrd/
        self.ligand_activity_thresholds = {'Kinase': 7.522879,
                                           'GPCR': 7,
                                           'Nuclear Receptor': 7,
                                           'Ion Channel': 5,
                                           'Non-IDG': 6}
        self.feasibility_decision_data_types = ['Affected pathway',
                                                'Animal model',
                                                # 'Genetic association',
                                                # 'Literature',
                                                # 'Rna expression',
                                                # 'Somatic mutation'
                                                ]
        self.colnms = []
        self.basicinfocols = []
        self.bucketcols = []
        self.gtexcols = []
        self.barresmousecols = []
        self.barreshumancols = []
        self.hpaenrichmentcols = []
        self.hpaexpressioncols = []
        self.riskfactorcols = []
        self.diseasecols = []
        self.jacksonlabcols = []
        self.cansarcols = []
        self.antitargetcols = []
        self.corefitnesscols = []
        self.druggabilitycols = []
        self.existingdrugscols = []
        self.pharoscols = []
        self.gwascols = []
        self.antibodyabilitycols = []
        self.ab_loc_cols = []
        self.literaturecols = []
        self.pathwayscols = []
        
        self.set_hardcoded_lists()
        # Set vars for some other data sources; we'll fill them on run (if they aren't set already)
        # These expression datasets HAVE to be filled in during startup, because (some of) their headers are used in the
        # setting up of disease profiles. Inefficient, but difficult to avoid without doing some very silly things. 
        self.bucketing_guide = None
        self.id_lookup_table = None
        self.locationdata = None
        self.barreslab_mouse = None
        self.bl_mouse_cell_types = []
        self.bl_mouse_scaler = StandardScaler()
        self.barreslab_human = None
        self.bl_human_cell_types = []
        self.bl_human_scaler = StandardScaler()
        self.gtexdata = None
        self.gtexdata_tissue_types = []
        self.gtex_scaler = StandardScaler()
        
        self.gwas_catalog = None
        self.ebi_drugebility_data = None
        self.ebi_drugebility_lookup = None
        self.core_fitness_data = None
        self.crispr_core_fitness_data = None
        self.ogee_essential_gene_data = None
        self.surfaceome_inclusion = None
        self.ecm_components = None
        self.gene_family_data = {}
        self.antitargets_data = None
        self.cansar_data = None
        self.HPO_genes_to_diseases = None
        self.HPO_diseases_to_genes = None
        self.HPO_annotations = {}
        self.HPO_annotations_scores = None
        self.serious_disease_phenotypes = []
        self.mouse_severe_phenotypes_patterns = []
        self.mgi_phenotypes_to_top_level_classifications = {}
        self.identity_scores = None
        self.evalues = None
        self.pathway_genes = {}
        self.pathway_xml = {}
        
        # Make lookup tables that can be dynamically filled as we go, reducing repeated querying. 
        self.pharos_ligand_data_library = {}
        self.drugebility_data_for_uniprot_ids = {}
        
        self.clinical_trial_data = {}
        self.opentargets_datatypes = {}
        self.pdb_to_uniprot = None
        self.pdb_data = None
        self.ppi_data = None
        self.hpa_expression_data = None
        self.drug_tox_data = None
        self.dgidb_interactions = None
        self.chembl_interactions = None
        self.iuphar_gene_ligand_interactions = None
        self.t3db_toxins = None
        self.t3db_moas = None
        self.t3db_drug_moas = None
        self.verseda_secretome = None
        self.literature = None
        self.withdrawn_withdrawn = None
        self.gnomad = None
        self.hpa_tox_lists = {}
        self.hpa_tox_tissues = ['Heart muscle', 'Liver', 'Kidney', 'Testis', 'Cervix, uterine', 'Ovary', 'Bone marrow',
                                'Lymph node']
        self.hpa_listclasses = []
        
        # Set dataframes that will be built up as output
        self.targets = None
        self.jackson_lab_data = None
        self.monogenic_disease_association_data = None
        self.polygenic_disease_association_data = None
        self.existing_drug_data = None
        self.gwas_data = None
        self.pathways = None
        
        # Set up connection objects
        self.mousemine = None
        self.server = None
        self.biomart_dataset = None
        self.drugebility_is_online = True
        self.ot = None
        
        # A few counters
        self.successful_requests = []
        self.druggable_count = 0
    
    def set_api_keys(self, api_keys):
        self.api_keys = api_keys
    
    def read_expression_datasets(self):
        # This slows down the startup of the program by reading in some datasets, but the hassle of dealing with the 
        # alternative is just not worth it. I need these headers to be available when picking tissues for disease
        # profiles. 
        
        # Barres lab data
        if not self.barreslab_mouse:
            self.status.emit("Reading Barres lab data (mouse)...")
            print("Reading Barres lab data (mouse)...")
            self.barreslab_mouse = read_excel("Data/barreslab_rnaseq (2).xlsx", sheet_name=None)
            self.bl_mouse_cell_types = self.barreslab_mouse['Raw Data'].columns.tolist()[2:]
            self.bl_mouse_scaler.fit(self.barreslab_mouse['Raw Data'][self.barreslab_mouse['Raw Data'].columns[2:]])
        
        if not self.barreslab_human:
            self.status.emit("Reading Barres lab data (human)...")
            print("Reading Barres lab data (human)...")
            self.barreslab_human = read_excel("Data/BarresLab_TableS4-HumanMouseMasterFPKMList.xlsx",
                                              sheet_name="Human data only",
                                              index_col=0)
            self.bl_human_scaler.fit(self.barreslab_human.values[2:])
            # Column names on that need a bit of work:
            cols = []
            hed = None
            for i, j in zip(self.barreslab_human.columns, self.barreslab_human.loc['Gene'].tolist()):
                if not match("Unnamed: [0-9]+", i):
                    hed = i
                cols.append(str(hed + ": " + j))
            self.barreslab_human.columns = cols
            self.bl_human_cell_types = list(self.barreslab_human.columns)
        
        # GTEx data
        try:
            if not self.gtexdata:
                self.status.emit("Reading GTEx data...")
                print("Reading GTEx data...")
                # self.gtexdata = read_table("Data/GTEx_Analysis_v6p_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct")
                self.gtexdata = read_csv("Data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.txt",
                                         low_memory=False, sep='\t')
                self.gtex_scaler.fit(self.gtexdata[self.gtexdata.columns[2:]])
                # GTEx data is organised by Ensembl ID, but it's got a version number on it. We're not interested in 
                # that right now. Get rid of it so we can do direct string equality comparisons. They're way easier to 
                # write. 
                # self.gtexdata['Ensembl ID'] = self.gtexdata['Name'].apply(lambda x: sub("\.[0-9]+", '', x))
                self.gtexdata['Ensembl ID'] = self.gtexdata['gene_id'].apply(lambda x: sub("\.[0-9]+", '', x))
                # Get a list of the tissue types listed in this file.
                self.gtexdata_tissue_types = [i for i in self.gtexdata.columns.tolist()
                                              if i not in ['Name', 'Ensembl ID', 'gene_id', 'Description']]
                self.colnms = list(set(self.colnms).union(set(self.gtex_tissues)))
                self.colnms = list(set(self.colnms).union(set(self.gtexdata_tissue_types)))
                # This is neat - to switch between HPA GTEx data and, uh, GTEx GTEx data, just swap between the columns 
                # requested here. Obviously, both sets can be included if deemed necessary. 
                # self.gtexcols.extend(gtex_tissues) 
                self.gtexcols.extend(self.gtexdata_tissue_types)
        except Exception as e:
            print("Exception in GTEX data read-in:")
            print(e)
            print_exc()
    
    def read_datasets(self):
        # A guide to bucketing (this just gets added to the output as a sheet)
        if not self.bucketing_guide:
            print("Reading bucketing guide...")
            self.bucketing_guide = read_excel("Data/bucket_guide.xlsx")
        
        # A fairly comprehensive basic ID lookup table, courtesy of GeneNames
        if not self.id_lookup_table:
            print("Reading ID lookup table...")
            self.id_lookup_table = read_csv("Data/genenames_all_gene_ids_lookup_table.tsv", sep='\t',
                                            converters={i: str for i in range(0, 100)})  # Read all cols as strings
            # Not bad, but that leaves us with empty strings, which we don't want. Replace them with NaN
            self.id_lookup_table = self.id_lookup_table.replace('', nan)
        
        # Subcellular location data
        if not self.locationdata:
            print("Reading subcellular locations...")
            self.status.emit("Reading subcellular location data...")
            self.locationdata = read_csv("Data/subcellular_location.csv")
        
        # GWAS catalog data
        if not self.gwas_catalog:
            print("Reading GWAS catalog...")
            self.status.emit("Reading GWAS Catalog data...")
            self.gwas_catalog = read_csv("Data/gwas_catalog_v1.0.1-associations_e88_r2017-04-24.tsv", sep='\t',
                                         low_memory=False)
        
        # EBI DrugEBIlity data
        if not self.ebi_drugebility_data:
            print("Reading DrugEBIlity...")
            self.ebi_drugebility_data = read_csv("Data/domain_drugebility.txt", sep='\t')
            self.ebi_drugebility_lookup = read_csv("Data/domain_details.txt", sep='\t')
        
        # List of core fitness genes
        if not self.core_fitness_data:
            print("Reading core fitness data...")
            self.status.emit("Reading core fitness genes data...")
            self.core_fitness_data = read_csv("Data/core_fitness_genes.csv")
        
        # List of CRISPR-screened core fitness genes
        if not self.crispr_core_fitness_data:
            print("Reading CRISPR screening core fitness data...")
            self.status.emit("Reading CRISPR-screened core fitness genes data...")
            self.crispr_core_fitness_data = read_csv("Data/Human_core_essential_gene_set_v2.txt", sep='\t',
                                                     names=["HGNC Name", "HGNC ID"])
        
        if not self.ogee_essential_gene_data:
            print("Reading OGEE...")
            self.status.emit("Reading OGEE essential genes data...")
            self.ogee_essential_gene_data = read_csv("Data/ogee_human_essential_genes.csv")
        
        # Surfaceome inclusion list
        if not self.surfaceome_inclusion:
            print("Reading surfaceome data...")
            self.status.emit("Reading surfaceome inclusion list data...")
            self.surfaceome_inclusion = read_excel("Data/Surfaceome_inclusion_lists.xlsx", sheet_name=None)
        
        # Extracellular matrix components list
        if not self.ecm_components:
            print("Reading extracellular matrix data...")
            self.status.emit("Reading extracellular matrix components list data...")
            self.ecm_components = read_excel("Data/corematrisome_hs.xls")
        
        # canSAR data (I downloaded everything I could get into a single file, in lieu of them providing an API)
        if not self.cansar_data:
            print("Reading CanSAR data...")
            self.status.emit("Reading canSAR data...")
            self.cansar_data = read_csv("Data/canSAR_all_human_genes_2017-12-11.csv")
            self.cansar_data.columns = [sub("_", " ", i) for i in list(self.cansar_data.columns)]
            extra_cansar_cols = [i for i in list(self.cansar_data.columns) if
                                 i not in ['search term', 'uniprot accession', 'gene name', 'description', 'synonyms',
                                           'location', 'keywords', 'cansar_link']]
            self.cansarcols = self.cansarcols + extra_cansar_cols
            for i in extra_cansar_cols:
                if i not in self.colnms:
                    self.colnms.append(i)
        
        # Read in disease phenotype-gene associations. These files don't actually need modifying from their original
        # format, so updating should be trivial. 
        if not self.HPO_genes_to_diseases:
            print("Reading HPO data...")
            self.status.emit("Reading Human Protein Ontology data...")
            self.HPO_genes_to_diseases = read_csv("Data/HPO_genes_to_diseases.txt", sep="\t", skiprows=1,
                                                  names=['EntrezID', 'HGNC Name', 'Disease'])
        
        if not self.HPO_diseases_to_genes:
            print("Reading HPO data (other way)...")
            self.status.emit("Reading Human Protein Ontology data...")
            self.HPO_diseases_to_genes = read_csv("Data/HPO_diseases_to_genes.txt", sep="\t",
                                                  skiprows=1, names=['Disease', 'EntrezID', 'HGNC Name'])
        
        if not self.pathway_genes:
            print("Reading WikiPathway data...")
            pathway_files = glob("Data/wikipathways-20181010-gpml-Homo_sapiens/*.gpml")
            for pfl in pathway_files:
                pathway_name = sub('.gpml', '', basename(pfl))
                # print(pathway_name)
                with open(pfl, encoding="utf8") as fd:
                    doc = xmltodict.parse(fd.read())
                labels = [i['@TextLabel'] for i in doc['Pathway']['DataNode'] if i.get('@Type') == "GeneProduct"]
                self.pathway_genes[pathway_name] = labels
                self.pathway_xml[pathway_name] = doc
        
        # Read a list of tissues we've picked out as being likely toxicity indicators from the Human Protein Atlas.
        # We can now auto-populate some lists that will help us keep this organised as desired.
        print("Reading toxicity data...")
        for fl in hpa_toxfiles:
            df = read_csv(fl, index_col=0)
            # dt = data type
            tis, dt = [sub('\+', ' ', i).capitalize() for i in split('[_\.]', basename(fl))[0:2]]
            self.hpa_listclasses.append(dt)
            if not self.hpa_tox_lists.get(tis):
                self.hpa_tox_lists[tis] = {}
            self.hpa_tox_lists[tis][dt] = df
        self.hpa_listclasses = sorted(uniq(self.hpa_listclasses))
        # Add tissue names to the Risk factors sheet 
        for i in self.hpa_tox_tissues:
            j = "HPA " + i
            if j not in self.riskfactorcols:
                self.riskfactorcols.append(j)
        # And the overall pandas frame columns! Super important. 
        for i in hpa_all_tissues:
            j = "HPA " + i
            if j not in self.colnms:
                self.colnms.append(j)
        
        if not self.HPO_annotations:
            print("Reading HPO annotations...")
            # This will involve a bit of somewhat dynamic programming. All good, as long as we do it carefully.
            for i in hpo_annotation_files:
                j = basename(i)
                j = sub("_diseases.csv", '', j)
                self.HPO_annotations[j] = read_csv(i, skiprows=1)
            # And add some column names
            for i in hpo_annotation_names:
                if i not in self.riskfactorcols:
                    self.riskfactorcols.append(i)
                if i not in self.colnms:
                    self.colnms.append(i)
        
        if not self.HPO_annotations_scores:
            print("Reading HPO annotaions scores...")
            self.HPO_annotations_scores = read_csv("Data/HPO_termname_lookup.csv")
        
        # Tag phenotypes with associated genes in the list of essential genes as being particularly serious.
        # We'll use this in safety classifications. 
        self.serious_disease_phenotypes = self.HPO_genes_to_diseases[self.HPO_genes_to_diseases['HGNC Name'].isin(
            self.core_fitness_data['HGNC Name'].tolist() + self.crispr_core_fitness_data['HGNC Name'].tolist()
        )]['Disease'].tolist()
        
        # Get identity scores (and e-values later, if necessary) for all genes in the human proteome.
        # I generated these using blastp. The JSON files get made in blast_handler.py
        print("Reading identity and e-value scores from all-vs-all BLASTp...")
        with open("Data/ident_matrix.json", 'r') as f:
            self.identity_scores = load(f)
        with open("Data/eval_matrix.json", 'r') as f:
            self.evalues = load(f)
        
        # Get PDB lookup table - links PDB codes to UniProt IDs.
        # I made this file myself, but it's based on https://www.uniprot.org/docs/pdbtosp.txt
        if not self.pdb_to_uniprot:
            print("Reading PDB data...")
            self.pdb_to_uniprot = read_csv("Data/pdb_to_uniprot.csv")
        
        # Drug-gene interaction database file full of interactions
        # There are at least three of these I can go after. DGIdb is one. ChEMBL is another. DrugBank is a third, though
        # that will have to be accessed differently. 
        if not self.dgidb_interactions:
            print("Reading DGIdb data...")
            self.dgidb_interactions = read_csv("Data/DGIdb_interactions.tsv", sep='\t')
        if not self.chembl_interactions:
            print("Reading ChEMBL interactions data...")
            self.chembl_interactions = read_csv("Data/chembl_drugtargets_named-18_13_46_02.csv",
                                                encoding="ISO-8859-1")
        
        # IUPHAR gene-ligand interaction database
        print("Reading IUPHAR data...")
        self.iuphar_gene_ligand_interactions = read_csv("Data/IUPHAR_ligand-gene_interactions.csv",
                                                        low_memory=False)
        self.iuphar_gene_ligand_interactions = self.iuphar_gene_ligand_interactions[
            self.iuphar_gene_ligand_interactions["target_species"] == "Human"]
        
        # Read these t3db things - toxicity database t3db.ca
        if not self.t3db_toxins:
            print("Reading T3DB toxicology data...")
            self.status.emit("Reading T3DB toxicology data...")
            self.t3db_toxins = read_json("Data/t3db_toxins.json")
        if not self.t3db_moas:
            self.t3db_moas = read_csv("Data/t3db_moas.csv")
        if not self.t3db_drug_moas:
            # Filter the methods of action down to only drugs
            drug_toxins = DataFrame(
                [row for i, row in self.t3db_toxins.iterrows() if 'Drug' in [j['type_name'] for j in row['types']]])
            self.t3db_drug_moas = self.t3db_moas[self.t3db_moas['Toxin T3DB ID'].isin(drug_toxins['title'])]
        
        # VerSeDa secretome genes - filtered by query to SecretomeP >= 0.9
        print("Reading VerSeDa data...")
        self.verseda_secretome = read_csv("Data/VerSeDa_HSapiens_secretedProteins.csv")
        
        # Read a list of withdrawn drugs from the WITHDRAWN database
        if not self.withdrawn_withdrawn:
            print("Reading WITHDRAWN database data...")
            self.withdrawn_withdrawn = read_csv("Data/withdrawn_withdrawn_compounds.csv")
        
        # Read Gnomad data
        if not self.gnomad:
            print("Reading Gnomad database data...")
            self.gnomad = read_csv("Data/gnomad.v2.1.1.lof_metrics.by_gene.txt", sep='\t')
        
        # We can now get Jackson lab mouse phenotype data, but it's organised in a form that is incompatible with the 
        # existing targets table. (Multiple lines per gene, basically - and no human reader-friendly means of 
        # reshaping it). Consequently, we have to create a new table, which we will deal with at the end. 
        self.jackson_lab_data = DataFrame(columns=self.jacksonlabcols)
        # OpenTargets (And OMIM too) likely return multiple associations per gene. Make a separate sheet, just like 
        # the Jackson lab data. 
        # In fact, make two sheets - one for rare/single-gene disease associations, another for common/polygenic
        # disease associations.
        self.monogenic_disease_association_data = DataFrame(columns=self.diseasecols)
        self.polygenic_disease_association_data = DataFrame(columns=self.diseasecols)
        # Also, make an 'existing drugs' table - that's extremely useful information. Also comes from OpenTargets. 
        self.existing_drug_data = DataFrame(columns=self.existingdrugscols)
        # Make an 'anitargets' sheet - targets with a similar protein sequence.
        self.antitargets_data = DataFrame(columns=self.antitargetcols)
        # GWAS data should be handled in a similar manner. 
        self.gwas_data = DataFrame(columns=self.gwascols)
        # PDB data also needs a sheet of this kind, but I'll be appending stuff to it, so it can be blank.
        self.pdb_data = DataFrame(columns=['series', 'GeneID', 'HGNC Name', 'Uniprot ID'])
        # Same again for protein-protein interaction data.
        self.ppi_data = DataFrame(columns=['series', 'GeneID', 'HGNC Name', 'Uniprot ID'])
        # And again for HPA RNA expression data
        self.hpa_expression_data = DataFrame(columns=self.hpaexpressioncols)
        # And again for drug toxicology data
        self.drug_tox_data = DataFrame(columns=['series', 'GeneID', 'HGNC Name'])
        # And again for literature list
        self.literature = DataFrame(columns=self.literaturecols)
        # Pathways needs a sheet...
        self.pathways = DataFrame(columns=self.pathwayscols)
        
        # Set these up as empty frames. They're our output. 
        # Working frame - we'll be adding stuff primarily to this.
        self.targets = self.input_data.copy()
        # Set empty columns up; it make everything easier later on. 
        for col in self.poss_existing_cols:
            if col not in self.targets.columns.tolist():
                self.targets[col] = None
        for i in self.colnms + self.hpa_predicted_secreted_protein_subclasses \
                + self.hpa_predicted_membrane_protein_subclasses:
            self.targets[i] = None
        # Strip whitespace characters from important columns (which may or may not currently contain data, but will at
        # least now consistently exist).
        for col in ['HGNC Name', 'GeneID', 'Uniprot ID']:
            self.targets[col] = self.targets[col].str.strip()
    
    def establish_connections(self):
        # Also, the connection to the MouseMine database needs to be set up.
        try:
            self.status.emit("Connecting to MouseMine...")
            self.mousemine = Service('www.mousemine.org/mousemine')
            # Set up a connection to BioMart too.
            self.status.emit("Connecting to BioMart...")
            self.server = BiomartServer("http://ensembl.org/biomart")
            # print("Checking BioMart connection...")
            # if server.is_alive:
            #     print("   Connected")
            # else:
            #     exit("   Unable to establish BioMart connection! Exiting.")
            self.biomart_dataset = self.server.datasets['hsapiens_gene_ensembl']
            self.status.emit("Connecting to OpenTargets...")
            
            # if not self.api_keys['OpenTargets_appName'] or len(str(self.api_keys['OpenTargets_appName'])) < 1:
            #     self.ot = OpenTargetsClient()
            # else:
            #     self.ot = OpenTargetsClient(auth_app_name=self.api_keys['OpenTargets_appName'],
            #                                 auth_secret=self.api_keys['OpenTargets_secret'])
            # Looks like they've removed this API key stuff entirely!
            self.ot = OpenTargetsClient()
            self.status.emit("Checking EBI DrugEBIlity status...")
            response = get("https://www.ebi.ac.uk/chembl/drugebility")
            if response.status_code != 200:
                self.warnings.emit("EBI Druggability query failed with code " + str(response.status_code))
                self.drugebility_is_online = False
                self.warnings.emit("EBI DrugEBIlity is not accessible.\nOther data will be downloaded as normal.")
            self.status.emit("Connections established")
        except Exception as e:
            print("Exception in connection establisher function:")
            print(e)
            print_exc()
            exit(1)
    
    def set_hardcoded_lists(self):
        # HPA protein classes.We'll use these shortly.
        self.toplevel = {'Enzymes', 'CD markers', 'Blood group antigen proteins', 'Nuclear receptors', 'Transporters',
                         'Ribosomal proteins', 'G-protein coupled receptors', 'Voltage-gated ion channels',
                         'Predicted membrane proteins', 'Predicted secreted proteins',
                         'Predicted intracellular proteins',
                         'Plasma proteins', 'Transcription factors', 'RNA polymerase related proteins',
                         'RAS pathway related proteins', 'Citric acid cycle related proteins'}
        self.secondlevel = {'ENZYME proteins', 'Peptidases', 'Kinases', 'Transporter channels and pores',
                            'Electrochemical Potential-driven transporters', 'Primary Active Transporters',
                            'Transport Electron Carriers', 'Accessory Factors Involved in Transport',
                            'Family 3 (C) receptors',
                            'GPCRs excl olfactory receptors', 'Adenosine and adenine nucleotide receptors',
                            'Chemokines and chemotactic factors receptors', 'Lysolipids receptors',
                            'Odorant/olfactory and gustatory receptors', 'Opsins', 'Serotonin receptors',
                            'Family 2 (B) receptors',
                            'Family T2R receptors (taste receptor GPCRs)', 'Family fz/smo receptors',
                            'Calcium-Activated Potassium Channels', 'CatSper and Two-Pore Channels',
                            'Cyclic Nucleotide-Regulated Channels', 'Inwardly Rectifying Potassium Channels',
                            'Transient Receptor Potential Channels', 'Two-P Potassium Channels',
                            'Voltage-Gated Calcium Channels',
                            'Voltage-Gated Potassium Channels', 'Voltage-Gated Sodium Channels',
                            'Prediction method-based',
                            '# TM segments-based', 'Secreted proteins predicted by MDSEC',
                            'SignalP predicted secreted proteins',
                            'Phobius predicted secreted proteins', 'SPOCTOPUS predicted secreted proteins',
                            'Yet undefined DNA-binding domains', 'Basic domains',
                            'Zinc-coordinating DNA-binding domains',
                            'Helix-turn-helix domains', 'Other all-alpha-helical DNA-binding domains',
                            'alpha-Helices exposed by beta-structures', 'Immunoglobulin fold',
                            'beta-Hairpin exposed by an alpha/beta-scaffold', 'beta-Sheet binding to DNA',
                            'beta-Barrel DNA-binding domains', }
        self.thirdlevel = {'Oxidoreductases', 'Transferases', 'Hydrolases', 'Lyases', 'Isomerase', 'Ligase',
                           'Aspartic-type peptidases', 'Cysteine-type peptidases', 'Metallopeptidases',
                           'Serine-type peptidases',
                           'Threonine-type peptidases', 'AGC Ser/Thr protein kinases', 'Tyr protein kinases',
                           'TKL Ser/Thr protein kinases', 'STE Ser/Thr protein kinases',
                           'RGC receptor guanylate cyclase kinases',
                           'NEK Ser/Thr protein kinases', 'CMGC Ser/Thr protein kinases', 'CK1 Ser/Thr protein kinases',
                           'CAMK Ser/Thr protein kinases', 'Atypical kinases', 'Membrane proteins predicted by MDM',
                           'MEMSAT3 predicted membrane proteins', 'MEMSAT-SVM predicted membrane proteins',
                           'Phobius predicted membrane proteins', 'SCAMPI predicted membrane proteins',
                           'SPOCTOPUS predicted membrane proteins', 'THUMBUP predicted membrane proteins',
                           'TMHMM predicted membrane proteins', 'GPCRHMM predicted membrane proteins',
                           '1TM proteins predicted by MDM', '9TM proteins predicted by MDM',
                           '8TM proteins predicted by MDM',
                           '7TM proteins predicted by MDM', '6TM proteins predicted by MDM',
                           '5TM proteins predicted by MDM',
                           '4TM proteins predicted by MDM', '3TM proteins predicted by MDM',
                           '2TM proteins predicted by MDM',
                           '>9TM proteins predicted by MDM'}
        self.hpa_predicted_membrane_protein_subclasses = ["Membrane proteins predicted by MDM",
                                                          "MEMSAT3 predicted membrane proteins",
                                                          "MEMSAT-SVM predicted membrane proteins",
                                                          "Phobius predicted membrane proteins",
                                                          "SCAMPI predicted membrane proteins",
                                                          "SPOCTOPUS predicted membrane proteins",
                                                          "THUMBUP predicted membrane proteins",
                                                          "TMHMM predicted membrane proteins",
                                                          "GPCRHMM predicted membrane proteins"]
        self.hpa_predicted_secreted_protein_subclasses = ["Secreted proteins predicted by MDSEC",
                                                          "SignalP predicted secreted proteins",
                                                          "Phobius predicted secreted proteins",
                                                          "SPOCTOPUS predicted secreted proteins"]
        
        self.poss_existing_cols = ['Uniprot ID', 'HGNC Name', 'GeneID', 'series', 'Approved Name', 'Previous Name',
                                   'Previous Symbols', 'Synonyms', 'Mouse Ensembl ID', 'Mouse Uniprot ID', 'MGI Symbol',
                                   'Mouse associated gene name', 'Entrez ID', 'Chromosome',
                                   'HCOP Rat', 'HCOP Worm', 'HCOP Mouse', 'HCOP Fly', 'HCOP Zebrafish', 'Gene OMIM ID']
        
        # This is not an ideal way of doing this, but it prevents some quite annoying load times and removes some 
        # annoying complexity. 
        self.gtex_tissues = ['Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)', 'Adrenal Gland',
                             'Artery - Aorta', 'Artery - Coronary', 'Artery - Tibial', 'Bladder', 'Brain - Amygdala',
                             'Brain - Anterior cingulate cortex (BA24)', 'Brain - Caudate (basal ganglia)',
                             'Brain - Cerebellar Hemisphere', 'Brain - Cerebellum', 'Brain - Cortex',
                             'Brain - Frontal Cortex (BA9)', 'Brain - Hippocampus', 'Brain - Hypothalamus',
                             'Brain - Nucleus accumbens (basal ganglia)', 'Brain - Putamen (basal ganglia)',
                             'Brain - Spinal cord (cervical c-1)', 'Brain - Substantia nigra',
                             'Breast - Mammary Tissue', 'Cells - EBV-transformed lymphocytes',
                             'Cells - Transformed fibroblasts', 'Cervix - Ectocervix', 'Cervix - Endocervix',
                             'Colon - Sigmoid', 'Colon - Transverse', 'Esophagus - Gastroesophageal Junction',
                             'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube',
                             'Heart - Atrial Appendage', 'Heart - Left Ventricle', 'Kidney - Cortex', 'Liver', 'Lung',
                             'Minor Salivary Gland', 'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary', 'Pancreas',
                             'Pituitary', 'Prostate', 'Skin - Not Sun Exposed (Suprapubic)',
                             'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum', 'Spleen', 'Stomach',
                             'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole Blood']
        self.barreslab_human_tissues = ['Human GBM / peri-tumor astrocytes', 'Human sclerotic hippocampi astrocytes',
                                        'Human fetal astrocytes', 'Human mature astrocytes', 'Human Neurons',
                                        'Human Oligodendrocytes', 'Human Microglia/Macrophage', 'Human Endothelial',
                                        'Human whole cortex']
        self.barreslab_mouse_tissues = ['cerebral cortex', 'Astrocytes', 'Neuron', 'Oligodendrocyte Precursor Cell',
                                        'Newly Formed Oligodendrocyte', 'Myelinating Oligodendrocytes', 'Microglia',
                                        'Endothelial cells']
        self.general_tissues = ['Adipose - Subcutaneous', 'Adipose - Visceral (Omentum)', 'Adrenal Gland',
                                'Artery - Aorta', 'Artery - Coronary', 'Artery - Tibial', 'Bladder', 'Brain - Amygdala',
                                'Brain - Anterior cingulate cortex (BA24)', 'Brain - Caudate (basal ganglia)',
                                'Brain - Cerebellar Hemisphere', 'Brain - Cerebellum', 'Brain - Cortex',
                                'Brain - Frontal Cortex (BA9)', 'Brain - Hippocampus', 'Brain - Hypothalamus',
                                'Brain - Nucleus accumbens (basal ganglia)', 'Brain - Putamen (basal ganglia)',
                                'Brain - Spinal cord (cervical c-1)', 'Brain - Substantia nigra',
                                'Breast - Mammary Tissue', 'Cells - EBV-transformed lymphocytes',
                                'Cells - Transformed fibroblasts', 'Cervix - Ectocervix', 'Cervix - Endocervix',
                                'Colon - Sigmoid', 'Colon - Transverse', 'Esophagus - Gastroesophageal Junction',
                                'Esophagus - Mucosa', 'Esophagus - Muscularis', 'Fallopian Tube',
                                'Heart - Atrial Appendage', 'Heart - Left Ventricle', 'Kidney - Cortex', 'Liver',
                                'Lung', 'Minor Salivary Gland', 'Muscle - Skeletal', 'Nerve - Tibial', 'Ovary',
                                'Pancreas', 'Pituitary', 'Prostate', 'Skin - Not Sun Exposed (Suprapubic)',
                                'Skin - Sun Exposed (Lower leg)', 'Small Intestine - Terminal Ileum', 'Spleen',
                                'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole Blood', 'Human Endothelial',
                                'Human whole cortex', 'cerebral cortex', ]
        self.cell_types = ['Human GBM / peri-tumor astrocytes', 'Human sclerotic hippocampi astrocytes',
                           'Human fetal astrocytes', 'Human mature astrocytes', 'Human Neurons',
                           'Human Oligodendrocytes', 'Human Microglia/Macrophage', 'Astrocytes', 'Neuron',
                           'Oligodendrocyte Precursor Cell', 'Newly Formed Oligodendrocyte',
                           'Myelinating Oligodendrocytes', 'Microglia', 'Endothelial cells']
        
        self.colnms = ['Input type', 'RNA class', 'Locus group', 'Feature type', 'Mammalian phenotype ID', 'MPID term',
                       'cerebral cortex', 'Astrocytes', 'Neuron', 'Oligodendrocyte Precursor Cell',
                       'Newly Formed Oligodendrocyte',
                       'Myelinating Oligodendrocytes', 'Microglia', 'Endothelial Cells',
                       # 'Consensus', 'LY druggability',
                       'EBI Tractable', 'EBI Druggable', 'EBI Ensemble',
                       'Domain ID', 'PDB', 'Druggable class', 'Pharos', 'ChEMBL drug', 'ChEMBL ligand',
                       'ChEMBL low-activity ligand', 'Is protein',
                       'Main location', 'Top level protein classes', 'Second level protein classes',
                       'Third level protein classes', 'Housekeeping gene',
                       'Has withdrawn drug', 'Withdrawn drug list',
                       'Membrane proteins predicted by MDM', 'GPCRHMM predicted membrane proteins', '# TM segments',
                       'Predicted secreted proteins', 'Mutational cancer driver genes',
                       'COSMIC somatic mutations in cancer genes', 'Gene family', 'Gene family ID', 'Ligand',
                       'Endogenous ligand', 'Rare diseases', 'ORPHANET', 'GWAS catalogue', 'PHEWAS catalogue',
                       'Core fitness gene', 'CRISPR-screened core fitness gene', 'OGEE human essential gene',
                       'HPA query', 'GTEX query', 'Location data query', 'DrugEBIlity query', 'Core fitness query',
                       'MGI query', 'Barres lab mouse query', 'Barres lab human query', 'Pharos query', 'Pharos ID',
                       'OpenTargets query', 'GWAS query', 'GPCR', 'Predicted membrane protein',
                       'Predicted secreted protein', 'VerSeDa secretome membership', 'Extracellular matrix component',
                       'Subcellular location summary', 'Subcellular location verification',
                       'Main subcellular locations', 'Additional subcellular locations', 'GOIDs',
                       'Surfaceome membership (human)', 'Surfaceome membership confidence (human)',
                       'Surfaceome membership (mouse)', 'Surfaceome membership confidence (mouse)',
                       'Functional summary', 'Phenotype risk score', 'Withdrawn drug ratio',
                       'Human GBM / peri-tumor astrocytes: 59yo surround',
                       'Human GBM / peri-tumor astrocytes: 59yo tumor core',
                       'Human GBM / peri-tumor astrocytes: 65yo tumor core',
                       'Human GBM / peri-tumor astrocytes: 64yo tumor core',
                       'Human sclerotic hippocampi astrocytes: 21yo',
                       'Human sclerotic hippocampi astrocytes: 22yo',
                       'Human sclerotic hippocampi astrocytes: 53yo - A',
                       'Human sclerotic hippocampi astrocytes: 53yo - B', 'Human fetal astrocytes: 18gw - 1',
                       'Human fetal astrocytes: 18gw - 2', 'Human fetal astrocytes: 18.5gw - 1',
                       'Human fetal astrocytes: 18.1gw', 'Human fetal astrocytes: 18.5gw - 2',
                       'Human fetal astrocytes: 18.5gw - 3', 'Human mature astrocytes: 8yo',
                       'Human mature astrocytes: 13yo', 'Human mature astrocytes: 16yo',
                       'Human mature astrocytes: 21yo', 'Human mature astrocytes: 22yo',
                       'Human mature astrocytes: 35yo', 'Human mature astrocytes: 47yo',
                       'Human mature astrocytes: 51yo', 'Human mature astrocytes: 53yo',
                       'Human mature astrocytes: 60yo', 'Human mature astrocytes: 63yo - 1',
                       'Human mature astrocytes: 63yo - 2', 'Human Neurons: 25yo', 'Human Oligodendrocytes: 22yoGC',
                       'Human Oligodendrocytes: 63yoGC - 1', 'Human Oligodendrocytes: 63yo GC - 2',
                       'Human Oligodendrocytes: 47yoO4', 'Human Oligodendrocytes: 63yoO4',
                       'Human Microglia/Macrophage: 45yo', 'Human Microglia/Macrophage: 51yo',
                       'Human Microglia/Macrophage: 63yo', 'Human Endothelial: 13yo', 'Human Endothelial: 47yo',
                       'Human whole cortex: 45yo', 'Human whole cortex: 63yo', 'Human whole cortex: 25yo',
                       'Human whole cortex: 53yo',
                       'Assays', 'Antibodies', 'Affected pathway', 'Animal model', 'Genetic association', 'Known drug',
                       'Literature', 'Pharos literature', 'Rna expression', 'Somatic mutation',
                       'Safety bucket', 'SM Druggability bucket', 'Feasibility bucket', 'AB-ability bucket',
                       'New modality bucket', 'Tissue engagement bucket', 'gnomAD pLI', 'gnomAD o/e',
                       'OpenTargets Ensemble', 'OpenTargets high quality compounds', 'Small molecule genome member'
                       ]
        self.basicinfocols = ['series', 'HGNC Name', 'GeneID', 'Approved Name', 'Previous Name', 'Previous Symbols',
                              'Synonyms', 'Uniprot ID', 'Entrez ID', 'Chromosome', 'Gene family', 'Gene family ID',
                              'Functional summary', 'Mouse Ensembl ID', 'Mouse Uniprot ID', 'MGI Symbol',
                              'HCOP Rat', 'HCOP Worm', 'HCOP Fly', 'HCOP Zebrafish', 'HCOP Mouse', 'Locus group',
                              'Is protein', 'RNA class', 'Top level protein classes']
        self.bucketcols = ['series', 'HGNC Name', 'GeneID', 'Safety bucket', 'SM Druggability bucket',
                           'Feasibility bucket', 'AB-ability bucket', 'New modality bucket', 'Tissue engagement bucket',
                           'Pharos']
        self.gtexcols = ['series', 'HGNC Name', 'GeneID', 'GTEX query']
        self.barresmousecols = ['series', 'HGNC Name', 'GeneID', 'Uniprot ID', 'Barres lab mouse query',
                                'cerebral cortex', 'Astrocytes', 'Neuron', 'Oligodendrocyte Precursor Cell',
                                'Newly Formed Oligodendrocyte', 'Myelinating Oligodendrocytes', 'Microglia',
                                'Endothelial Cells']
        self.barreshumancols = ['series', 'HGNC Name', 'GeneID', 'Uniprot ID', 'Barres lab human query',
                                'Human GBM / peri-tumor astrocytes: 59yo surround',
                                'Human GBM / peri-tumor astrocytes: 59yo tumor core',
                                'Human GBM / peri-tumor astrocytes: 65yo tumor core',
                                'Human GBM / peri-tumor astrocytes: 64yo tumor core',
                                'Human sclerotic hippocampi astrocytes: 21yo',
                                'Human sclerotic hippocampi astrocytes: 22yo',
                                'Human sclerotic hippocampi astrocytes: 53yo - A',
                                'Human sclerotic hippocampi astrocytes: 53yo - B', 'Human fetal astrocytes: 18gw - 1',
                                'Human fetal astrocytes: 18gw - 2', 'Human fetal astrocytes: 18.5gw - 1',
                                'Human fetal astrocytes: 18.1gw', 'Human fetal astrocytes: 18.5gw - 2',
                                'Human fetal astrocytes: 18.5gw - 3', 'Human mature astrocytes: 8yo',
                                'Human mature astrocytes: 13yo', 'Human mature astrocytes: 16yo',
                                'Human mature astrocytes: 21yo', 'Human mature astrocytes: 22yo',
                                'Human mature astrocytes: 35yo', 'Human mature astrocytes: 47yo',
                                'Human mature astrocytes: 51yo', 'Human mature astrocytes: 53yo',
                                'Human mature astrocytes: 60yo', 'Human mature astrocytes: 63yo - 1',
                                'Human mature astrocytes: 63yo - 2', 'Human Neurons: 25yo',
                                'Human Oligodendrocytes: 22yoGC', 'Human Oligodendrocytes: 63yoGC - 1',
                                'Human Oligodendrocytes: 63yo GC - 2', 'Human Oligodendrocytes: 47yoO4',
                                'Human Oligodendrocytes: 63yoO4', 'Human Microglia/Macrophage: 45yo',
                                'Human Microglia/Macrophage: 51yo', 'Human Microglia/Macrophage: 63yo',
                                'Human Endothelial: 13yo', 'Human Endothelial: 47yo', 'Human whole cortex: 45yo',
                                'Human whole cortex: 63yo', 'Human whole cortex: 25yo', 'Human whole cortex: 53yo']
        self.hpaenrichmentcols = ['series', 'HGNC Name', 'GeneID', 'HPA query']
        self.hpaexpressioncols = ['series', 'HGNC Name', 'GeneID', 'HPA query', 'Adipose tissue', 'Adrenal gland', 
                                  'Amygdala', 'Appendix', 'B-cells', 'Basal ganglia', 'Bone marrow', 'Breast', 
                                  'Cerebellum', 'Cerebral cortex', 'Cervix, uterine', 'Colon', 'Corpus callosum', 
                                  'Dendritic cells', 'Ductus deferens', 'Duodenum', 'Endometrium 1', 'Epididymis', 
                                  'Esophagus', 'Fallopian tube', 'Gallbladder', 'Granulocytes', 'Heart muscle', 
                                  'Hippocampal formation', 'Hypothalamus', 'Kidney', 'Liver', 'Lung', 'Lymph node', 
                                  'Midbrain', 'Monocytes', 'NK-cells', 'Olfactory region', 'Ovary', 'Pancreas', 
                                  'Parathyroid gland', 'Pituitary gland', 'Placenta', 'Pons and medulla', 'Prostate', 
                                  'Rectum', 'Retina', 'Salivary gland', 'Seminal vesicle', 'Skeletal muscle', 'Skin 1', 
                                  'Small intestine', 'Smooth muscle', 'Spinal cord', 'Spleen', 'Stomach 1', 'T-cells', 
                                  'Testis', 'Thalamus', 'Thymus', 'Thyroid gland', 'Tongue', 'Tonsil', 'Total PBMC', 
                                  'Urinary bladder', 'Vagina']
        self.riskfactorcols = ['series', 'HGNC Name', 'GeneID',
                               'Phenotype risk score', 'Withdrawn drug ratio',
                               'Mutational cancer driver genes', 'COSMIC somatic mutations in cancer genes',
                               'Housekeeping gene', 'gnomAD pLI', 'gnomAD o/e']
        self.pathwayscols = ['series', 'HGNC Name', 'GeneID', 'Uniprot ID', 'Pathway name', 'Pathway description']
        self.diseasecols = ['series', 'HGNC Name', 'GeneID', 'Uniprot ID', 'OpenTargets query', 'Approved Name',
                            'Source', 'Disease name', 'Therapeutic areas', 'Disease ID', 'Association score',
                            'Genetic association', 'Direct association']
        self.jacksonlabcols = ['series', 'HGNC Name', 'GeneID', 'Approved name', 'Input', 'Input Type',
                               'MGI Gene/Marker ID', 'Mouse associated gene symbol', 'Feature type', 'MP ID', 'Term',
                               'Top-level phenotype', 'Description']
        self.cansarcols = ['series', 'HGNC Name', 'GeneID']
        self.antitargetcols = ['series', 'HGNC Name', 'GeneID', 'Uniprot ID', 'Antitarget gene name',
                               'Antitarget Uniprot ID', 'Protein seq. identity (%)']
        self.corefitnesscols = ['series', 'HGNC Name', 'GeneID', 'Core fitness query', 'Core fitness gene',
                                'CRISPR-screened core fitness gene', 'OGEE human essential gene']
        self.druggabilitycols = ['series', 'HGNC Name', 'GeneID', 'DrugEBIlity query',
                                 # 'Consensus', 'LY druggability',
                                 'EBI Tractable', 'EBI Druggable', 'EBI Ensemble', 'OpenTargets Ensemble',
                                 # 'Top Druggable Domain Ensemble',
                                 'Domain ID', 'PDB', 'Is protein', 'Small molecule genome member',
                                 # 'Druggable class', 'Pharos', 'ChEMBL drug', 'ChEMBL ligand',
                                 'Ligand', 'Endogenous ligand', 'Main location', 'Top level protein classes',
                                 'Second level protein classes',
                                 'Third level protein classes', 'Membrane proteins predicted by MDM',
                                 'GPCRHMM predicted membrane proteins', '# TM segments',
                                 'OpenTargets high quality compounds',
                                 # 'Predicted secreted proteins'
                                 'Has withdrawn drug', 'Withdrawn drug list'
                                 ]
        self.existingdrugscols = ['series', 'HGNC Name', 'GeneID', 'Disease', 'Association score', 'Drug name',
                                  'Molecule type', 'Action type', 'Mechanism of action', 'Associated', 'Resource score',
                                  'Max clinical phase', 'Boxed warning', 'Drug ID', 'Withdrawn',
                                  'Withdrawn reason', 'Withdrawn country', 'Withdrawn year']
        self.pharoscols = ['series', 'HGNC Name', 'GeneID', 'Pharos query', 'Pharos ID', 'Druggable class', 'Pharos',
                           'ChEMBL drug', 'ChEMBL ligand', 'ChEMBL low-activity ligand', 'Pharos literature']
        self.gwascols = ['series', 'HGNC Name', 'GeneID', 'Source', 'Study name', 'ID', 'Phenotypes', 'Highest P-value',
                         'Num. markers', 'link']
        self.antibodyabilitycols = ['series', 'HGNC Name', 'GeneID', 'Uniprot ID', 'GPCR', 'Predicted membrane protein',
                                    'Predicted secreted protein', 'VerSeDa secretome membership',
                                    'Extracellular matrix component',
                                    'Subcellular location summary', 'Subcellular location verification',
                                    'Main subcellular locations', 'Additional subcellular locations', 'GOIDs',
                                    'Surfaceome membership (human)', 'Surfaceome membership confidence (human)',
                                    'Surfaceome membership (mouse)',
                                    'Surfaceome membership confidence (mouse)', ] + \
                                   self.hpa_predicted_secreted_protein_subclasses + \
                                   self.hpa_predicted_membrane_protein_subclasses
        self.feasibilitycols = ['series', 'HGNC Name', 'GeneID', 'Uniprot ID', "Assays", "Antibodies",
                                # 'Affected pathway', 'Animal model',
                                # 'Genetic association', 
                                # 'Known drug', 
                                'Literature',
                                # 'Rna expression', 'Somatic mutation', 
                                ] + self.feasibility_decision_data_types
        self.literaturecols = ['series', 'GeneID', 'HGNC Name', 'Uniprot ID', 'Authors', 'DOI', 'Date', 'Journal ref',
                               'Title']
        self.mouse_severe_phenotypes_patterns = ["premature death", "decreased survivor rate", "lethality"]
    
    def run(self):
        try:
            if self.input_data is None:
                # We can't proceed without a list of targets supplied.
                print("No input data. It's all gone terribly wrong somewhere.")
                exit(1)
            
            # Read in in-house/pre-downloaded datasets
            # These have been separated into two functions because some are necessary in some cases (needed for setting 
            # up disease profiles), so must occasionally be loaded before the others. 
            self.read_expression_datasets()
            self.read_datasets()
            # Set up connections to external databases where necessary
            self.establish_connections()
            
            print("Starting")
            success_count = 0
            self.druggable_count = 0
            num_targets = len(self.targets)
            for ind, target in self.targets.iterrows():
                # A bit of housekeeping first: we may or may not have useful things like gene name, Ensembl ID, 
                # and UniProt ID. (Presumably we'll have at least one of those; what are we even working with if not?)
                # If any are missing, they should be filled in. MAY need to shuffle this around based on preference for 
                # getting one type of info from another. 
                target = self.get_basic_annotations(ind, target)
                
                # Bear in mind that we may not have an HGNC name. (This will be the case for RNAs, clones, etc.)
                # If we don't, it needs to be handled in a robust manner. 
                displayname = self.set_display_name(target)
                # self.status.emit("Downloading data for target " + displayname + "...")
                self.emit_download_status(target, displayname, " ")
                
                # Has an "Approved name" been assigned? If not, use genenames.org to find it. 
                # Look for synonyms while we're at it...
                print("genenames")
                self.get_genenames_annotations(ind, target, displayname)
                
                # The NCBI has some useful annotations too. Go and get 'em. 
                print("ncbi")
                self.get_ncbi_annotations(ind, target, displayname)
                self.get_ncbi_assays(ind, target, displayname)
                
                # Get names of homologs in other model species using BioMart.
                # (We can get TONS more from BioMart if we want to, but this is all we really need right now).
                print("biomart")
                self.get_biomart_annotations(ind, target, displayname)
                
                # Check for interactions with ligands in IUPHAR data
                print("iuphar")
                self.get_ligands(ind, target, displayname)
                
                # Human Protein Atlas.
                print("hpa")
                self.get_hpa_data(ind, target, displayname)
                
                # A minor complication here: There appear to be two, slightly different sets of GTEx data out in the 
                # wild. One is accessible via the HPA, the other directly via the GTEx portal. The above code dealt with 
                # the former, the code below with the latter (which I read in earlier). Eventually we'll pick one or the 
                # other, but for now, we'll prepare both.
                print("gtex")
                self.get_gtex_data(ind, target, displayname)
                
                print("hpa enrichment")
                self.get_hpa_enrichment(ind, target, displayname)
                
                # Look up location data (and level of confidence) from locationdata
                print("location")
                self.get_location_data(ind, target, displayname)
                
                # gnomAD
                print("gnomad")
                self.get_gnomad_data(ind, target, displayname)
                
                # Druggability. 
                print("drugebility")
                # if self.drugebility_is_online:
                #     self.get_druggability_data(ind, target, displayname)
                self.get_ebi_drugebility(ind, target, displayname)
                
                # canSAR
                print("canSAR")
                self.get_cansar_data(ind, target, displayname)
                
                # PDB
                print("pdb")
                self.get_pdb_structures(ind, target, displayname)
                
                # Anti-targets
                print("anti-targets")
                self.get_antitargets(ind, target, displayname)
                
                # WikiPathways
                print("pathways")
                self.get_pathways_data(ind, target, displayname)
                
                # Is this gene in the core fitness gene list?
                print("core fitness")
                self.get_core_fitness_data(ind, target, displayname)
                
                # Is this gene in the extracellular matrix gene list?
                print("ecm")
                self.get_extracellular_matrix_data(ind, target, displayname)
                
                # Is this gene in the VerSeDa secretome gene list?
                print("verseda secretome")
                self.get_verseda_secretome_data(ind, target, displayname)
                
                # Jackson lab mouse data.
                # There is likely to be more than one record per target here, so we'll need to use a different 
                # arrangement of output sheet. 
                print("jaxlab")
                self.get_jackson_lab_data(ind, target, displayname)
                
                # Barres lab data.
                # This is basically just merging matching rows from any columns in that dataset with matching names.
                # Also, if there is no mouse gene symbol attributed, we can't do anything here. 
                print("barres lab")
                self.get_barres_lab_data(ind, target, displayname)
                
                # Pharos
                print("pharos")
                self.get_pharos_data(ind, target, displayname)
                
                # Disease association data 
                # This will be a bit tricky because, as with the mouse lab data, there is likely
                # to be more than one association per gene. May need to create another sheet again. 
                # We get a lot of other stuff in here too - literature associated a disease to a target, and also some
                # information about drugs.
                print("opentargets")
                self.get_disease_association_data(ind, target, displayname)
                
                # Risk factors - some stuff that relates to safety, but which isn't so easy to directly make a bucketing
                # call based on. This has to happen after we've used OpenTargets.
                print("risk factors")
                self.get_risk_factors(ind, target, displayname)
                
                # Drug safety data
                # A follow-on from OpenTargets, where we got a list of existing drugs.
                # We're specifically checking for the existence of boxed warnings. 
                print("fda")
                self.get_drug_safety_data(ind, target, displayname)
                
                # Get GWAS results from GWAS Central and GWAS Catalog
                print("gwas")
                self.get_gwas_data(ind, target, displayname)
                
                # Get AB-ability data (or, at least, the parts that haven't been sorted out elsewhere yet)
                print("ab-ability")
                self.get_antibodyability(ind, target, displayname)
                
                # Get protein-protein interaction data
                # This adds a HUGE amount of stuff to the output. Disable for now.
                print("protein-protein")
                self.get_protein_protein_interactions(ind, target, displayname)
                
                # Get drug toxicology data
                print("drug tox")
                self.get_drug_toxicology(ind, target, displayname)
                
                # Finish it off by updating the progress bar
                if self.successful_requests:
                    success_count += 1
                self.progbar_update.emit(int(((ind + 1) / num_targets) * 100))
                print("done with that target!\n-----")
            
            # This where to think about doing stuff that comes after downloading data is essentially 
            # complete. That is, basic analysis and the like. I.e., bucketing!
            self.bucketing()
            
            self.status.emit("Writing output file...")
            print("TOTAL:\n   " + str(success_count) + " genes found in DrugEBIlity, of " + str(len(self.targets)) +
                  "\n   " + str(self.druggable_count) + " of those have druggable domains")
            
            sheets_by_name = self.gather_output_sheets()
            
            # print("Got " + str(len(sheets_by_name)) + " sheets")
            # for i in sheets_by_name:
            #     print(i)
            #     print(sheets_by_name[i])
            
            # return sheets_by_name
            
            print("Returning data!")
            self.got_data.emit(sheets_by_name)
        except Exception as e:
            print("Exception in base miner function:")
            print(e)
            print_exc()
            print("This is broad enough that I can't auto-resolve it, unfortunately. ")
            exit(1)
    
    def emit_download_status(self, target, displayname, text):
        # A convenience function for writing a status message to the GUI.
        self.status.emit(
            "Downloading data for target " + str(target['series']) + ": " + displayname + "... (" + text + ")")
    
    def emit_bucket_status(self, target):
        # A convenience function for writing a status message to the GUI. (A different one to the download status).
        self.status.emit(
            "Assigning buckets for target " + str(target['series']))
    
    def get_basic_annotations(self, ind, target):
        # Take care of filling in HGNC Name, Ensembl ID (GeneID) and UniProt ID in the most effective order.
        if not isnull(target['HGNC Name']):
            # If we have a gene name...
            print("Starting with HGNC name " + target['HGNC Name'])
            if isnull(target['GeneID']):
                # ... but no HGNC Name, get Ensembl ID
                ensembl_id = self.get_ensembl_id_for_gene_name(target['HGNC Name'])
                if ensembl_id:
                    self.targets.at[ind, 'GeneID'] = ensembl_id
                    target['GeneID'] = ensembl_id  # Make it available inside the loop too
                    print("\t-> " + ensembl_id)
            if isnull(target['Uniprot ID']):
                # ... but no Uniprot ID, get Uniprot ID
                uniprot_id = self.get_uniprot_id_for_gene_name(target['HGNC Name'])
                if uniprot_id:
                    self.targets.at[ind, 'Uniprot ID'] = uniprot_id
                    target['Uniprot ID'] = uniprot_id  # Make it available inside the loop too
                    print("\t-> " + uniprot_id)
        elif not isnull(target['GeneID']):
            # Otherwise, if we have an Ensembl ID...
            print("Starting with Ensembl ID " + target['GeneID'])
            if isnull(target['HGNC Name']):
                # ... but no name, get a name
                gene_name = self.get_gene_name_for_ensembl_id(target['GeneID'])
                if gene_name:
                    self.targets.at[ind, 'HGNC Name'] = gene_name
                    target['HGNC Name'] = gene_name  # Make it available inside the loop too
                    print("\t-> " + gene_name)
                else:
                    print("\t-> No gene name found")
            if isnull(target['Uniprot ID']):
                # ... but no Uniprot ID, get Uniprot ID
                uniprot_id = self.get_uniprot_id_for_ensembl_id(target['GeneID'])
                if uniprot_id:
                    self.targets.at[ind, 'Uniprot ID'] = uniprot_id
                    target['Uniprot ID'] = uniprot_id  # Make it available inside the loop too
                    print("\t-> " + uniprot_id)
                else:
                    print("\t-> no uniprot ID found")
        elif not isnull(target['Uniprot ID']):
            # Otherwise, if we have a Uniprot ID...
            print("Starting with Uniprot ID " + target['Uniprot ID'])
            if isnull(target['HGNC Name']):
                # ... but no name, get a name
                gene_name = self.get_gene_name_for_uniprot_id(target['Uniprot ID'])
                if gene_name:
                    self.targets.at[ind, 'GeneID'] = gene_name
                    target['HGNC Name'] = gene_name  # Make it available inside the loop too
                    print("\t-> " + gene_name)
            if isnull(target['GeneID']):
                # ... but no GeneID, get Ensembl ID 
                # ensembl_id = get_ensembl_id_for_uniprot_id(target['Uniprot ID']) 
                # This one is tricky; I haven't found a way of directly converting one to the other via 
                # APIs and such; only via manual forms at most. Instead, we'll have to look up the Ensembl ID 
                # from the name (which is OK, because we just tried to fill it in). 
                if not isnull(target['HGNC Name']):
                    ensembl_id = self.get_ensembl_id_for_gene_name(target['HGNC Name'])
                    if ensembl_id:
                        self.targets.at[ind, 'GeneID'] = ensembl_id
                        target['GeneID'] = ensembl_id  # Make it available inside the loop too
                        print("\t-> " + ensembl_id)
        return target
    
    def set_display_name(self, target):
        # A systematic way of choosing which ID of a target to display in the GUI.
        if target['HGNC Name'] not in [None, nan, 'nan']:
            return target['HGNC Name']
        elif target['GeneID'] not in [None, nan, 'nan']:
            return "[unnamed target " + target['GeneID'] + "]"
        elif target['Uniprot ID'] not in [None, nan, 'nan']:
            return "[unnamed target " + target['Uniprot ID'] + "]"
        else:
            exit("Something has gone very wrong with gene input.")
            # print(name)
    
    def get_genenames_annotations(self, ind, target, displayname):
        # Get some annotations - family membership and the like - that's directly accessible via basic identifiers.
        try:
            self.emit_download_status(target, displayname, "GeneNames annotations")
            # Now that I get all this stuff from a file rather than web queries, I can quite radically simplify this:
            df = None
            if target['HGNC Name'] not in [None, nan, 'nan']:
                df = self.id_lookup_table[self.id_lookup_table['Approved Symbol'] == target['HGNC Name']]
            elif target['GeneID'] not in [None, nan, 'nan']:
                df = self.id_lookup_table[self.id_lookup_table['Ensembl Gene ID'] == target['GeneID']]
            elif target['Uniprot ID'] not in [None, nan, 'nan']:
                df = self.id_lookup_table[
                    self.id_lookup_table['UniProt ID(supplied by UniProt)'] == target['Uniprot ID']]
            else:
                exit("Extremely incorrect input data found; should have been caught long before this!")
            
            if len(df):
                # The names of things in the lookup table and my table often don't quite match. Sort that out here.
                #                    GeneNames                       Mine
                pairs_of_things = [["Approved Name", "Approved Name"],
                                   ["Entrez Gene ID", "Entrez ID"],
                                   ["Previous Symbols", "Previous Symbols"],
                                   ["Previous Name", "Previous Name"],
                                   ["Synonyms", "Synonyms"],
                                   ["OMIM ID(supplied by OMIM)", "Gene OMIM ID"],
                                   ["Gene Family Name", "Gene family"],
                                   ["Gene Family ID", "Gene family ID"],
                                   ["Chromosome", "Chromosome"],
                                   ["Locus Group", "Locus group"]]
                # Now that that's taken care of, we can just go and pull those out and slot them into self.targets
                for i in pairs_of_things:
                    gn_name = i[0]
                    my_name = i[1]
                    value = df[gn_name].tolist()[0]
                    if isnull(value):
                        value = None
                    self.targets.at[ind, my_name] = value
                    target[my_name] = value
                    # A bit of extra stuff I nearly forgot about: gene family data. (I'll keep querying this for now).
                    if my_name == "Gene family ID" and value:
                        for gfid in self.split_gene_families(value):
                            self.get_protein_family_info(gfid)
        except Exception as e:
            print("Exception in GeneNames function:")
            print(e)
            print_exc()
            exit()
            # print("Retrying in 10 seconds...")
            # sleep(10)
            # self.get_genenames_annotations(ind, target, displayname)
    
    @staticmethod
    def split_gene_families(gene_family_id):
        # This bit of code gets used several times, so it should be a function.
        # Why do we need it? Because there may be multiple family IDs in this string!
        return gene_family_id.split("|")
    
    def get_protein_family_info(self, gene_family_id):
        # We may repeatedly query the same family here. In the interests of both speed and being nice, let's not 
        # actually do that. 
        # IMPORTANT: Some genes can belong to more than one family. But I've dealt with that elsewhere.
        url = None
        from pandas.errors import ParserError
        try:
            if gene_family_id not in self.gene_family_data:
                url = "https://www.genenames.org/cgi-bin/genefamilies/set/" + gene_family_id + "/download/" + \
                      self.gene_family_query_level
                self.gene_family_data[gene_family_id] = read_csv(url, sep='\t')
        except ParserError as pe:
            print("Parser error; get to the bottom of this.")
            print(pe)
            print(url)
            print(gene_family_id)
            print(type(gene_family_id))
            print_exc()
            exit()
        except Exception as e:
            print("Exception in GeneNames gene family function:")
            print(e)
            print(url)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            self.get_protein_family_info(gene_family_id)
    
    def get_biomart_annotations(self, ind, target, displayname):
        # Query BioMart for more advanced IDs, including homologs and some RNA stuff.
        try:
            if target['GeneID'] not in [None, nan, 'nan']:
                # self.status.emit("Downloading data for target " + displayname + "... (BioMart annotations)")
                self.emit_download_status(target, displayname, "BioMart annotations")
                attrs = ['celegans_homolog_associated_gene_name',
                         'celegans_homolog_ensembl_gene',
                         'drerio_homolog_associated_gene_name',
                         'drerio_homolog_ensembl_gene',
                         'rnorvegicus_homolog_associated_gene_name',
                         'rnorvegicus_homolog_ensembl_gene',
                         'dmelanogaster_homolog_associated_gene_name',
                         'dmelanogaster_homolog_ensembl_gene',
                         'mmusculus_homolog_associated_gene_name',
                         'mmusculus_homolog_ensembl_gene']
                response = self.biomart_dataset.search({'filters': {'ensembl_gene_id': target['GeneID']},
                                                        'attributes': attrs})
                biomrt = read_csv(StringIO("\n".join([i.decode('utf-8') for i in response.iter_lines()])), sep='\t',
                                  names=attrs)
                # Great; that now needs to be written into appropriate columns. 
                # Make these lists unique. We seem to get duplicates quite often.
                flylist = [str(i) for i in uniq(biomrt['dmelanogaster_homolog_associated_gene_name'].dropna().tolist())]
                if isnull(target['HCOP Fly']) and flylist:
                    self.targets.at[ind, 'HCOP Fly'] = ", ".join(flylist)
                ratlist = [str(i) for i in uniq(biomrt['rnorvegicus_homolog_associated_gene_name'].dropna().tolist())]
                if isnull(target['HCOP Rat']) and ratlist:
                    self.targets.at[ind, 'HCOP Rat'] = ", ".join(ratlist)
                wormlist = [str(i) for i in uniq(biomrt['celegans_homolog_associated_gene_name'].dropna().tolist())]
                if isnull(target['HCOP Worm']) and wormlist:
                    self.targets.at[ind, 'HCOP Worm'] = ", ".join(wormlist)
                fishlist = [str(i) for i in uniq(biomrt['drerio_homolog_associated_gene_name'].dropna().tolist())]
                if isnull(target['HCOP Zebrafish']) and fishlist:
                    self.targets.at[ind, 'HCOP Zebrafish'] = ", ".join(fishlist)
                mouselist = [str(i) for i in uniq(biomrt['mmusculus_homolog_ensembl_gene'].dropna().tolist())]
                if isnull(target['Mouse Ensembl ID']) and mouselist:
                    self.targets.at[ind, 'Mouse Ensembl ID'] = ", ".join(mouselist)
                    # Do this to prevent getting rate-limited by Ensembl when doing multiple queries sequentially
                    mouse_uniprot_list = []
                    for enid in mouselist:
                        mouse_uniprot_id = self.get_uniprot_id_for_ensembl_id(enid)
                        if mouse_uniprot_id:
                            mouse_uniprot_list.append(mouse_uniprot_id)
                        sleep(1)
                    # mouse_uniprot_list = list(
                    #     filter(None, [self.get_uniprot_id_for_ensembl_id(i, wait=True) for i in mouselist]))
                    self.targets.at[ind, 'Mouse Uniprot ID'] = ", ".join(mouse_uniprot_list)
                    target['Mouse Uniprot ID'] = mouse_uniprot_list
                mouselist = [str(i) for i in uniq(biomrt['mmusculus_homolog_associated_gene_name'].dropna().tolist())]
                if isnull(target['HCOP Mouse']) and mouselist:
                    self.targets.at[ind, 'HCOP Mouse'] = ", ".join(mouselist)
                    target['HCOP Mouse'] = ", ".join(mouselist)  # Because we'll actually be using this one later.
                print("Mouse homolog names: " + ", ".join(mouselist))
            
            # RNA classes information
            # I'm not 100% sure on this, but it's the closest I've been able to find.
            # It's got to be its own query; BioMart INSISTS on it.
            print("rna classes")
            if target['GeneID'] not in [None, nan, 'nan']:
                attrs = ['transcript_biotype']
                response = self.biomart_dataset.search({'filters': {'ensembl_gene_id': target['GeneID']},
                                                        'attributes': attrs})
                rna = read_csv(StringIO("\n".join([i.decode('utf-8') for i in response.iter_lines()])), sep='\t',
                               names=attrs)
                if isnull(target['RNA class']) and rna['transcript_biotype'].tolist():
                    self.targets.at[ind, 'RNA class'] = ", ".join(rna['transcript_biotype'].tolist())
            
            # GO terms?
            print("GO terms")
            print("    (work in progress)")
        except Exception as e:
            print("Exception in BioMart function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            self.get_biomart_annotations(ind, target, displayname)
    
    def get_ncbi_annotations(self, ind, target, displayname):
        # Get NCBI IDs
        try:
            if target['Entrez ID'] not in [None, nan, 'nan']:
                self.emit_download_status(target, displayname, "NCBI annotations")
                url = "https://www.ncbi.nlm.nih.gov/gene/" + str(target['Entrez ID']) + "?report=xml&format=text"
                response = get(url)
                data = xmltodict.parse(response.text)
                if response.status_code not in [200, 404]:
                    self.warnings.emit("NCBI query failed with code " + str(response.status_code))
                else:
                    if data:
                        if 'pre' in data:
                            # For some reason, we may have to do this twice
                            data = xmltodict.parse(data['pre'])
                        if 'Entrezgene' in data:
                            if 'Entrezgene_summary' in data['Entrezgene']:
                                self.targets.at[ind, 'Functional summary'] = data['Entrezgene']['Entrezgene_summary']
                                target['Functional summary'] = data['Entrezgene']['Entrezgene_summary']
        except Exception as e:
            print("Exception in NCBI annotations function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            self.get_ncbi_annotations(ind, target, displayname)
    
    def get_ncbi_assays(self, ind, target, displayname):
        # Get some information about bioassays for activity on a given gene. Counts towards feasibility.
        try:
            if target['HGNC Name']:
                # self.status.emit("Downloading data for target " + displayname + "... (NCBI assays)")
                self.emit_download_status(target, displayname, "NCBI assays")
                url = "https://www.ncbi.nlm.nih.gov/pcassay/?term=" + str(target['HGNC Name']) + "&format=text"
                print(url)
                response = get(url)
                if response.status_code == 200:
                    # The NCBI assay search tries to do something clever here, but it actually makes things a bit more 
                    # difficult for us. When there is only one result from a search, it just redirects us straight to 
                    # the relevant PubChem link. When this happens, the returned text is explicitly html, not xml.
                    if search("DOCTYPE html", response.text, IGNORECASE):
                        self.targets.at[ind, 'Assays'] = 1
                    else:
                        data = xmltodict.parse(response.text)
                        if data:
                            if data.get('pre'):
                                n = len(data['pre'].split("\n\n"))
                                if n == 20:
                                    n = "20+"
                                self.targets.at[ind, 'Assays'] = n
        
        except Exception as e:
            print("Exception in NCBI assays function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            self.get_ncbi_assays(ind, target, displayname)
    
    def get_ligands(self, ind, target, displayname):
        self.emit_download_status(target, displayname, "IUPHAR ligands")
        df = DataFrame(columns=['endogenous'])  # We've got to have this column, at least
        if target['HGNC Name']:
            df = self.iuphar_gene_ligand_interactions[
                self.iuphar_gene_ligand_interactions['target_gene_symbol'] == target['HGNC Name']]
        elif target['Uniprot ID']:
            df = self.iuphar_gene_ligand_interactions[
                self.iuphar_gene_ligand_interactions['target_uniprot'] == target['Uniprot ID']]
        self.targets.at[ind, 'Ligand'] = len(df)
        self.targets.at[ind, 'Endogenous ligand'] = len(df[df['endogenous'] == 't'])
    
    def get_hpa_data(self, ind, target, displayname):
        try:
            if target['GeneID'] not in [None, nan, 'nan']:
                self.emit_download_status(target, displayname, "Human Protein Atlas")
                # Get xml from the human protein atlas concerning this gene
                url = "http://www.proteinatlas.org/" + target['GeneID'] + ".xml"
                print("\t" + url)
                self.targets.at[ind, 'HPA query'] = target['GeneID']
                target['HPA query'] = target['GeneID']
                # This try/except is here because there can be cases where a gene isn't in the HPA. 
                data = None
                response = get(url)
                if response.status_code not in [200, 404]:
                    self.warnings.emit("HPA query failed with code " + str(response.status_code))
                else:
                    try:
                        data = xmltodict.parse(response.text)
                    except ExpatError as e:
                        print("  XML retrieval failed for this gene")
                    if data:
                        # pprint(data['proteinAtlas']['entry']['proteinClasses'])
                        if isinstance(data['proteinAtlas']['entry']['proteinClasses']['proteinClass'], list):
                            classes = set(
                                [i['@name'] for i in data['proteinAtlas']['entry']['proteinClasses']['proteinClass']])
                        else:
                            classes = set(data['proteinAtlas']['entry']['proteinClasses']['proteinClass']['@name'])
                        # Top-level protein classes (list of values)
                        if classes.intersection(self.toplevel):
                            self.targets.at[ind, 'Top level protein classes'] = ", ".join(
                                classes.intersection(self.toplevel))
                        # Second-level protein classes (list of values)
                        if classes.intersection(self.secondlevel):
                            self.targets.at[ind, 'Second level protein classes'] = ", ".join(
                                classes.intersection(self.secondlevel))
                        # Third-level protein classes (list of values)
                        if classes.intersection(self.thirdlevel):
                            self.targets.at[ind, 'Third level protein classes'] = ", ".join(
                                classes.intersection(self.thirdlevel))
                        # # Membrane proteins predicted by MDM (true/false)
                        # if 'Membrane proteins predicted by MDM' in classes:
                        #     self.targets.set_value(index=ind, col='Membrane proteins predicted by MDM', value=True)
                        # # GPCRHMM predicted membrane proteins (true/false)
                        # if 'GPCRHMM predicted membrane proteins' in classes:
                        #     self.targets.set_value(index=ind, col='GPCRHMM predicted membrane proteins', value=True)
                        # I now do this by collecting a full list of HPA's membrane protein/secreted protein prediction 
                        # methods - it comes in useful when doing AB-ability later.
                        for mpc in [i for i in classes if i in self.hpa_predicted_membrane_protein_subclasses]:
                            self.targets.at[ind, mpc] = True
                        for spc in [i for i in classes if i in self.hpa_predicted_secreted_protein_subclasses]:
                            self.targets.at[ind, spc] = True
                        # # TM segments (number of segments, based on max class)
                        tm_segs = [i for i in classes if search("TM proteins predicted by MDM", i)]
                        if tm_segs:
                            if '>9TM proteins predicted by MDM' in tm_segs:
                                self.targets.at[ind, '# TM segments'] = ">9"
                            else:
                                # If the >9 class ident is not present, we want the first character from the first class 
                                # ident when sorted in reverse order. That's the number of TM domains. 
                                self.targets.at[ind, '# TM segments'] = sorted(tm_segs, reverse=True)[0][0]
                        # Predicted secreted proteins (true/false)
                        if 'Predicted secreted proteins' in classes:
                            self.targets.at[ind, 'Predicted secreted proteins'] = True
                        # Mutational cancer driver genes (true/false)
                        if 'Mutational cancer driver genes' in classes:
                            self.targets.at[ind, 'Mutational cancer driver genes'] = True
                        # COSMIC somatic mutations  in cancer genes (true/false)
                        if 'COSMIC somatic mutations in cancer genes' in classes:
                            self.targets.at[ind, 'COSMIC somatic mutations in cancer genes'] = True
                        # Is protein - is in SwissProt (Yes + evidence level/No)
                        for i in data['proteinAtlas']['entry']['proteinEvidence']['evidence']:
                            if i['@source'] == "UniProt":
                                self.targets.at[ind, 'Is protein'] = "UniProt - " + i['@evidence']
                        # Get GTEX expression data for this protein   
                        # SWITCHING OFF FOR NOW - HPA have rearranged their returned data significantly.
                        # gtex = data['proteinAtlas']['entry']['rnaExpression']['data']
                        # List tissues present in this data:
                        # [i['tissue'] for i in gtex if 'tissue' in i]
                        # SWITCHING OFF FOR NOW - HPA have rearranged their returned data significantly.
                        # for i in gtex:
                        #     if 'tissue' in i:
                        #         tis = i['tissue']
                        #         if tis in self.gtex_tissues:
                        #             self.targets.at[ind, tis] = i['level']['@tpm']
                        # We can also get some data relating to localisation from the HPA.
                        # Get HPA expression data for this gene.
                        # This is different from, though similar to, the GTEX data I've commented out above. 
                        hpa_exp_vals = {j.get('tissue')['#text']: [k.get('@expRNA') for k in j.get('level') if
                                                                   k.get('@type') == "RNAExpression"]
                                        for j in [i.get('data') for i in
                                                  data.get('proteinAtlas').get('entry').get('rnaExpression')][0]}
                        for i in ['series', 'HGNC Name', 'GeneID', 'HPA query']:
                            hpa_exp_vals[i] = target[i]
                        # That may need some further tweaking to make it more robust. I'll try it and find out shortly.
                        self.hpa_expression_data = concat([self.hpa_expression_data, 
                                                           DataFrame.from_dict(hpa_exp_vals)], 
                                                          sort=False,
                                                          ignore_index=True)
                        
                        # I'll need to do some stuff that's a bit tricky here, but it's perfectly OK as long as I'm 
                        # careful about it.
                        # Some columns will be static, but others will be generated as new data comes in. They'll all be
                        # listed somewhere accessible though. It would also be a good idea to sort them alphabetically. 
                        # Easy stuff to start with:
                        # Is this a GPCR?
                        if 'G-protein coupled receptors' in classes:
                            self.targets.at[ind, 'GPCR'] = True
                        # Is this a predicted membrane protein (any method)?
                        if 'Predicted membrane proteins' in classes:
                            self.targets.at[ind, 'Predicted membrane protein'] = True
                        # Is this a predicted secreted protein (any method)?
                        if 'Predicted secreted proteins' in classes:
                            self.targets.at[ind, 'Predicted secreted protein'] = True
                        
                        # Is there any further localisation data available to us?
                        # GO term IDs are available here. Strip them out and hold on to them; we might use them later. 
                        if 'cellExpression' in data['proteinAtlas']['entry']:
                            cell_localisation_data = data['proteinAtlas']['entry']['cellExpression']
                            # actual values = cell_localisation_data['data']['location'][possible array]['#text']
                            loc_summary = cell_localisation_data['summary']
                            loc_verification = cell_localisation_data['verification']['#text']
                            # MAY have to deal with lists in the verification. Keep an eye on it. 
                            loc_main = []
                            loc_additional = []
                            goids = []
                            if isinstance(cell_localisation_data['data']['location'], list):
                                loc_main = [i['#text'] for i in cell_localisation_data['data']['location']
                                            if i['@status'] == "main"]
                                loc_additional = [i['#text'] for i in cell_localisation_data['data']['location']
                                                  if i['@status'] != "main"]
                                goids = [i['@GOId'] for i in cell_localisation_data['data']['location']
                                         if i['@status'] != "main"]
                            else:
                                loc_main.append(cell_localisation_data['data']['location']['#text'])
                                goids.append(cell_localisation_data['data']['location']['@GOId'])
                            self.targets.at[ind, 'Subcellular location summary'] = loc_summary
                            self.targets.at[ind, 'Subcellular location verification'] = loc_verification
                            self.targets.at[ind, 'Main subcellular locations'] = ", ".join(loc_main)
                            self.targets.at[ind, 'Additional subcellular locations'] = ", ".join(loc_additional)
                            self.targets.at[ind, 'GOIDs'] = ", ".join(goids)
                            # These may be useful later on; keep them handy.
                            target['Main subcellular locations'] = loc_main
                            target['Additional subcellular locations'] = loc_additional
                            target['GOIDs'] = goids
                            
                            # SWITCHING OFF FOR NOW - HPA have rearranged their returned data significantly.
                            # # Is this a housekeeping gene? (I.e., is it expressed in all tissues?)
                            # if data.get('proteinAtlas').get('entry').get('rnaExpression').get('rnaTissueCategory').get(
                            #         '#text') == "Expressed in all":
                            #     self.targets.at[ind, 'Housekeeping gene'] = True
                            
                            # We can also get a bit of information about assays and antibodies.
                            # This looks handy, but I don't see how to use it yet. 
                            # if 'antibody' in data['proteinAtlas']['entry']:
                            #     assay_availability = data['proteinAtlas']['entry']['antibody']
        except Exception as e:
            print("Exception in HPA function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            self.get_hpa_data(ind, target, displayname)
    
    def get_gtex_data(self, ind, target, displayname):
        if target['GeneID']:
            # self.status.emit("Downloading data for target " + displayname + "... (GTEX)")
            self.emit_download_status(target, displayname, "GTEX")
            self.targets.at[ind, 'GTEX query'] = target['GeneID']
            if target['GeneID'] in self.gtexdata['Ensembl ID'].tolist():
                gtex_hit = self.gtexdata.loc[self.gtexdata['Ensembl ID'] == target['GeneID']]
                for col in self.gtexdata_tissue_types:
                    self.targets.at[ind, col] = gtex_hit[col].values[0]
    
    def get_location_data(self, ind, target, displayname):
        # Subcellular location, that is.
        if target['HGNC Name'] or target['GeneID']:
            srch = 'GeneID'  # Field to search within targets
            mtch = 'Gene'  # Field in subcellular location data to search
            if target['HGNC Name']:
                srch = 'HGNC Name'
                mtch = 'Gene name'
            # self.status.emit("Downloading data for target " + displayname + "... (Subcellular localisation)")
            self.emit_download_status(target, displayname, "Subcellular localisation")
            self.targets.at[ind, 'Location data query'] = target[srch]
            if target[srch] in self.locationdata['Gene name'].unique():
                loc_conf = self.locationdata.loc[self.locationdata[mtch] == target[srch]]['Reliability'].values[0]
                location = self.locationdata.loc[self.locationdata[mtch] == target[srch]][loc_conf].values[0]
                self.targets.at[ind, 'Main location'] = str(location) + " (" + str(loc_conf) + ")"
    
    def get_ebi_drugebility(self, ind, target, displayname):
        # A replacement for the previous hacked-together web scraper thing for EBI's DrugEBIlity. 
        # This is now fairly straightforward because, thankfully, they have released the underlying data!
        # Unfortunately they've done that because the service is shutting down.
        try:
            self.emit_download_status(target, displayname, "DrugEBIlity")
            self.targets.at[ind, 'DrugEBIlity query'] = target['Uniprot ID']
            tuid = target['Uniprot ID']
            druggable = self.query_ebi_drugebility_data(tuid)
            if len(druggable) > 0:
                returnable_data = {'DrugEBIlity query': tuid,
                                   'EBI Tractable': druggable['TRACTABLE'].iloc[0],
                                   'EBI Druggable': druggable['DRUGGABLE'].iloc[0],
                                   'EBI Ensemble': druggable['ENSEMBLE'].iloc[0],
                                   'Domain ID': druggable['PX_NUMBER'].iloc[0],
                                   'PDB': druggable['PDB_CODE'].iloc[0]}
                self.drugebility_data_for_uniprot_ids[tuid] = returnable_data
                if returnable_data['EBI Druggable'] > 0.8:
                    self.druggable_count += 1
                for field in returnable_data:
                    self.targets.at[ind, field] = returnable_data[field]
        except Exception as e:
            print("Exception in DrugEBIlity function:")
            print(e)
            print_exc()
    
    def query_ebi_drugebility_data(self, tuid):
        pxs = self.ebi_drugebility_lookup[self.ebi_drugebility_lookup['ACCESSION'] == tuid]['PX_NUMBER'].tolist()
        df = self.ebi_drugebility_data[self.ebi_drugebility_data['PX_NUMBER'].isin(pxs)]
        # druggable = df[df['DRUGGABLE'] == 1]
        druggable = df.sort_values(by=["DRUGGABLE", "TRACTABLE", "ENSEMBLE"], ascending=False)
        return druggable
    
    def get_druggability_data(self, ind, target, displayname):
        self.emit_download_status(target, displayname, "DrugEBIlity")
        self.targets.at[ind, 'DrugEBIlity query'] = target['Uniprot ID']
        # target_uniprot_ids = reviewed[reviewed['Gene name'] == name]['Entry'] 
        self.successful_requests = []
        # Remember, we may get more than one of these per gene in the future...
        target_uniprot_ids = []
        if target['Uniprot ID']:
            target_uniprot_ids.append(target['Uniprot ID'])
        for tuid in target_uniprot_ids:
            drugebility = self.query_drugebility(tuid)
            if drugebility:
                for field in drugebility:
                    self.targets.at[ind, field] = drugebility[field]
    
    def query_drugebility(self, tuid):  # tuid = target uniprot id
        if self.drugebility_data_for_uniprot_ids.get(tuid):
            return self.drugebility_data_for_uniprot_ids[tuid]
        try:
            # self.targets.set_value(index=ind, col='DrugEBIlity query', value=tuid)
            # This is currently being troublesome, so let's try to not do that. 
            session = Session()
            retry = Retry(connect=5, backoff_factor=0.5)
            adapter = HTTPAdapter(max_retries=retry)
            session.mount('http://', adapter)
            session.mount('https://', adapter)
            response = session.get("https://www.ebi.ac.uk/chembl/drugebility/protein/structural/" + tuid)
            
            if response.status_code != 200:
                # self.warnings.emit("EBI Druggability query failed with code " + str(response.status_code))
                print("EBI Druggability query failed with code " + str(response.status_code))
                if response.status_code != 404:
                    sleep(10)
                    print("Re-attempting...")
                    self.query_drugebility(tuid)
            else:
                html = response.content
                # Can we pattern-match the phrase "Ave. Druggable" anywhere in that html? If so, we have the data 
                # we're looking for. 
                if search("Ave. Druggable", str(html)):
                    self.successful_requests.append(html)
                    soup = BeautifulSoup(html, "lxml")
                    tables = soup.find_all('table', {"class": "chembl_btable"})
                    # That picks out 2 tables. We want tables[1] right now.
                    # This cuts out the sub-table we used to use (also in tables[1], coincidentally), by telling pandas 
                    # to ignore the first 6 lines. We then use the next line as a header. 
                    # We also get some empty columns on the end, because of html page structure. Use .dropna 
                    # Also, sort descending by Ensemble score.
                    df = read_html(str(tables[0]), skiprows=6, header=0)[0] \
                        .dropna(axis=1) \
                        .sort_values(by='Ensemble', ascending=False)
                    
                    # Restrict the list to druggable domains only. (There may be none).
                    druggable = df[df['Druggable'] == 1]
                    if len(druggable) > 0:
                        self.druggable_count += 1
                        # This is where things start getting interesting. What should we do if we find multiple 
                        # druggable domains tied for the highest Ensemble score? Currently we just pick the first. 
                        # print("  Tractable = " + str(df[df[0] == 'Ave. Tractable'][1].values[0]))
                        # self.targets.set_value(index=ind, col='EBI Tractable', value=druggable['Tractable'].iloc[0])
                        # self.targets.set_value(index=ind, col='EBI Druggable', value=druggable['Druggable'].iloc[0])
                        # self.targets.set_value(index=ind, col='EBI Ensemble', value=druggable['Ensemble'].iloc[0])
                        # # self.targets.set_value(index=ind, col='Top Druggable Domain Ensemble',
                        # #                        value=druggable['Ensemble'].iloc[0])
                        # self.targets.set_value(index=ind, col='Domain ID', value=druggable['Domain ID'].iloc[0])
                        # self.targets.set_value(index=ind, col='PDB', value=druggable['PDB'].iloc[0])
                        returnable_data = {'DrugEBIlity query': tuid,
                                           'EBI Tractable': druggable['Tractable'].iloc[0],
                                           'EBI Druggable': druggable['Druggable'].iloc[0],
                                           'EBI Ensemble': druggable['Ensemble'].iloc[0],
                                           'Domain ID': druggable['Domain ID'].iloc[0],
                                           'PDB': druggable['PDB'].iloc[0]}
                        self.drugebility_data_for_uniprot_ids[tuid] = returnable_data
                        return returnable_data
        except Exception as e:
            print("Exception in DrugEBIlity function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            self.query_drugebility(tuid)
    
    def get_cansar_data(self, ind, target, displayname):
        # We can't query canSAR directly, but we can query a file I made from it. 
        try:
            if target['HGNC Name']:
                self.emit_download_status(target, displayname, "canSAR")
                cansar = self.cansar_data[self.cansar_data['gene name'] == target['HGNC Name']]
                # I think if I just drop a few columns containing info I already have, then we're already pretty close 
                # to having a usable sheet here. 
                # I set those columns up when I read in data. 
                if len(cansar) > 0:
                    for col in [i for i in self.cansarcols if i not in ['series', 'HGNC Name', 'GeneID']]:
                        self.targets.at[ind, col] = cansar[col].values[0]
        except Exception as e:
            print("Exception in canSAR function:")
            print(e)
            print_exc()
            exit("NOPEING OUT")
    
    def get_pathways_data(self, ind, target, displayname):
        try:
            if target['HGNC Name']:
                self.emit_download_status(target, displayname, "WikiPathways")
                pw = DataFrame(columns=self.pathwayscols)
                pwlist = [i for i in self.pathway_genes if target['HGNC Name'] in self.pathway_genes[i]]
                # pwlist = [i for i in pwlist if i in self.pathway_xml]
                pw['Pathway name'] = pwlist
                # pathway_desc = []
                # for xml in [self.pathway_xml[j] for j in pwlist]:
                #     if xml.get('Pathway'):
                #         if xml['Pathway'].get('Biopax'):
                #             if xml['Pathway']['Biopax'].get('bp:PublicationXref'):
                #                 
                #     pathway_desc.append(xml['Pathway']['Biopax']['bp:PublicationXref']['bp:TITLE']['#text'])
                
                # pathway_desc = [
                #     i['Pathway']['Biopax']['bp:PublicationXref']['bp:TITLE']['#text'] if i.get('Pathway').get(
                #         'Biopax').get('bp:PublicationXref').get('bp:TITLE').get('#text') else None for i in
                #         [self.pathway_xml[j] for j in pwlist]]
                
                xmlpath = ['Pathway', 'Biopax', 'bp:PublicationXref', 'bp:TITLE', '#text']
                pathway_desc = [deep_get(i, xmlpath) for i in [self.pathway_xml[j] for j in pwlist]]
                
                pw['Pathway description'] = pathway_desc
                for i in self.pathwayscols[:4]:
                    pw[i] = target[i]
                self.pathways = self.pathways.append(pw, ignore_index=True)
        except Exception as e:
            print("Exception in pathways function:")
            print(e)
            print_exc()
            exit()
    
    def get_core_fitness_data(self, ind, target, displayname):
        # Simple - check if a gene falls into any of these lists of core/essential genes. 
        if target['HGNC Name']:
            self.emit_download_status(target, displayname, "Core fitness genes")
            self.targets.at[ind, 'Core fitness query'] = target['HGNC Name']
            print("\tcore fitness")
            if target['HGNC Name'] in self.core_fitness_data['HGNC Name'].tolist():
                self.targets.at[ind, 'Core fitness gene'] = True
            print("\tCRISPR-screened core fitness")
            if target['HGNC Name'] in self.crispr_core_fitness_data['HGNC Name'].tolist():
                self.targets.at[ind, 'CRISPR-screened core fitness gene'] = True
            
            if target['HGNC Name'] in self.ogee_essential_gene_data['symbols'].tolist():
                val = self.ogee_essential_gene_data[self.ogee_essential_gene_data['symbols'] == target['HGNC Name']][
                    'essentiality consensus'].values[0]
                self.targets.at[ind, 'OGEE human essential gene'] = val
            elif target['GeneID']:
                if target['GeneID'] in self.ogee_essential_gene_data['locus'].tolist():
                    val = self.ogee_essential_gene_data[self.ogee_essential_gene_data['locus'] == target['GeneID']][
                        'essentiality consensus'].values[0]
                    self.targets.at[ind, 'OGEE human essential gene'] = val
        elif target['GeneID']:
            if target['GeneID'] in self.ogee_essential_gene_data['locus'].tolist():
                val = self.ogee_essential_gene_data[self.ogee_essential_gene_data['locus'] == target['GeneID']][
                    'essentiality consensus'].values[0]
                self.targets.at[ind, 'OGEE human essential gene'] = val
    
    def get_hpa_enrichment(self, ind, target, displayname):
        try:
            if target['HGNC Name']:
                print("hpa enrichment")
                self.emit_download_status(target, displayname, "HPA Enrichment")
                for t in hpa_all_tissues:
                    print("\t" + t)
                    outcomes = []
                    for enrichment in self.hpa_listclasses:
                        if target['HGNC Name'] in self.hpa_tox_lists[t][enrichment]['Gene'].tolist():
                            print("\t\t" + enrichment + " (length " + str(
                                len(self.hpa_tox_lists[t][enrichment]['Gene'].tolist())) + ")")
                            outcomes.append(enrichment)
                    if outcomes:
                        self.targets.at[ind, "HPA " + t] = ", ".join(outcomes)
        except Exception as e:
            print("Exception in HPA enrichment function:")
            print(e)
            print_exc()
            exit()
    
    def get_risk_factors(self, ind, target, displayname):
        # Check if our target has phenotypes that fall into any of the lists of phenotypes from the HPO corresponding to 
        # diseases affecting certain regions of the body, and calculate a score based upon that.
        try:
            if target['HGNC Name']:
                self.emit_download_status(target, displayname, "Risk factors")
                diseases = self.HPO_genes_to_diseases[self.HPO_genes_to_diseases['HGNC Name'] == target['HGNC Name']]
                score = 0
                for i in hpo_annotation_names:
                    df = self.HPO_annotations[i]
                    if list(set(diseases['Disease'].tolist()).intersection(df['Disease id'])):
                        self.targets.at[ind, i] = True
                        score += self.HPO_annotations_scores[self.HPO_annotations_scores['Term name'] == i][
                            'Withdrawn drug count'].values[0]
                self.targets.at[ind, 'Phenotype risk score'] = score
                # We have another column to add here as well. Specifically, we want to get the ratio of withdrawn drugs
                # to accepted drugs.
                drugs = self.existing_drug_data[self.existing_drug_data['HGNC Name'] == target['HGNC Name']]
                withdrawn_drugs = drugs[drugs['Withdrawn'] == True]
                ratio = 0
                if len(drugs):  # Avoid division by 0
                    ratio = float(len(withdrawn_drugs)) / float(len(drugs))
                self.targets.at[ind, 'Withdrawn drug ratio'] = ratio
                # This is a similar, but simpler, check. Cross-reference against a list of most strongly-associated 
                # genes for the most tox-prone tissues from the HPA.
                # for tissue in self.hpa_tox_tissues:
                #     outcomes = []
                #     for enrichment in self.hpa_listclasses:
                #         if target['HGNC Name'] in self.hpa_tox_lists[tissue][enrichment]['Gene'].tolist():
                #             outcomes.append(enrichment)
                #     if outcomes:
                #         self.targets.set_value(index=ind, col="HPA " + tissue, value=", ".join(outcomes))
                # Now dealt with elsewhere (self.get_hpa_enrichment())
        except Exception as e:
            print("Exception in risk factors function:")
            print(e)
            print_exc()
            exit()
    
    def get_extracellular_matrix_data(self, ind, target, displayname):
        # Another easy one - just look the gene up in the extracellular matrix data.
        if target['HGNC Name']:
            self.emit_download_status(target, displayname, "Extracellular matrix genes")
            print("\textracellular matrix")
            if target['HGNC Name'] in self.ecm_components['Gene Symbol'].tolist():
                self.targets.at[ind, 'Extracellular matrix component'] = True
    
    def get_verseda_secretome_data(self, ind, target, displayname):
        if target['GeneID']:
            self.emit_download_status(target, displayname, "VerSeDa secretome genes")
            if target['GeneID'] in self.verseda_secretome['Gene'].tolist():
                print("\tsecreted")
                self.targets.at[ind, 'VerSeDa secretome membership'] = True
        elif target['Uniprot ID']:
            self.emit_download_status(target, displayname, "VerSeDa secretome genes")
            if target['Uniprot ID'] in self.verseda_secretome['Uniprot'].tolist():
                print("\tsecreted")
                self.targets.at[ind, 'VerSeDa secretome membership'] = True
    
    def get_jackson_lab_data(self, ind, target, displayname):
        # Query MGI, the Jackson Lab's remarkably comprehensive database system.
        # This is best done with mouse homologs, if we have any. 
        try:
            if target['HCOP Mouse']:
                print("Jaxlab query: " + target['HCOP Mouse'])
                self.targets.at[ind, 'MGI query'] = target['HGNC Name']
                self.emit_download_status(target, displayname, "Jackson lab/MGI")
                template = self.mousemine.get_template('_Feature_Phenotype')
                template2 = self.mousemine.get_template('Term_Ancestors')  # We'll use this later
                rows = template.rows(B={"op": "LOOKUP", "value": target['HCOP Mouse']})
                # template = self.mousemine.get_template('HGene_MPhenotype')  # This doesn't quite do what we need.
                # rows = template.rows(A={"op": "LOOKUP", "value": target['HGNC Name']})
                buffer = DataFrame(columns=self.jacksonlabcols)
                for row in rows:
                    # row["subject.symbol"] is the mouse gene symbol. Record it - it's useful.
                    self.targets.at[ind, 'MGI Symbol'] = row["subject.symbol"]
                    target['MGI Symbol'] = row["subject.symbol"]
                    new_mouse_row = {'series': target['series'],
                                     'GeneID': target['GeneID'],
                                     'HGNC Name': target['HGNC Name'],
                                     'Approved Name': target['Approved Name'],
                                     "Input": target['HCOP Mouse'],
                                     "Input Type": "current symbol",
                                     "MGI Gene/Marker ID": row["subject.primaryIdentifier"],
                                     "Mouse associated gene symbol": row["subject.symbol"],
                                     "Feature type": row["subject.sequenceOntologyTerm.name"],
                                     "MP ID": row["ontologyTerm.identifier"],
                                     "Term": row["ontologyTerm.name"],
                                     "Description": row['evidence.comments.description']}
                    # We can add a bit more than that, if we do another query.
                    # Specifically, we can add top-level hierarchical classes.
                    # I believe this can raise unaccounted-for errors, though, so we need to handle that. 
                    # May be that I'm submitting too many repeated queries to MGI. Try getting only novel stuff.
                    print(row["ontologyTerm.identifier"])
                    
                    if row["ontologyTerm.identifier"] in self.mgi_phenotypes_to_top_level_classifications:
                        new_mouse_row['Top-level phenotype'] = self.mgi_phenotypes_to_top_level_classifications[
                            row["ontologyTerm.identifier"]]
                    else:
                        # sleep(1)  # JaxLab clearly imposes a query limit. Let's try to stay under it. 
                        
                        
                        # Alright, I tried being nice, but that didn't cut it. Time to stop caring.
                        def do_this_query():
                            try:
                                return template2.rows(A={"op": "LOOKUP", "value": row["ontologyTerm.identifier"]})
                            except Exception:
                                print("BEAST MODE")
                                # sleep(1)
                                do_this_query()
                        
                        phenrows = do_this_query()
                        
                        # phenrows = template2.rows(A={"op": "LOOKUP", "value": row["ontologyTerm.identifier"]})
                        parents = [i['OntologyTerm.parents.name'] for i in phenrows]
                        cunning_merge_thing = [i for i in high_level_phenotypes if
                                               [j for j in parents if match(i, j)]]
                        if cunning_merge_thing:
                            # This lets us handle rare cases where a phenotype is in multiple top-level classes
                            new_mouse_row['Top-level phenotype'] = ", ".join(cunning_merge_thing)
                            self.mgi_phenotypes_to_top_level_classifications[
                                row["ontologyTerm.identifier"]] = new_mouse_row['Top-level phenotype']
                    buffer = buffer.append(new_mouse_row, ignore_index=True)
                
                if len(buffer.index) == 0:
                    new_mouse_row = {'series': target['series'],
                                     'GeneID': target['GeneID'],
                                     'HGNC Name': target['HGNC Name'],
                                     'Approved Name': target['Approved Name'],
                                     "Input": target['HCOP Mouse'],
                                     "Input Type": "current symbol",
                                     "MGI Gene/Marker ID": "No phenotypes found"}
                    self.jackson_lab_data = self.jackson_lab_data.append(new_mouse_row, ignore_index=True)
                else:
                    # Do this all in one go at the end, in order to avoid inserting duplicate rows if we handle an error
                    # halfway through the list of rows...
                    for i, row in buffer.iterrows():
                        self.jackson_lab_data = self.jackson_lab_data.append(row, ignore_index=True)
            else:
                new_mouse_row = {'series': target['series'],
                                 'GeneID': target['GeneID'],
                                 'HGNC Name': target['HGNC Name'],
                                 'Approved Name': target['Approved Name'],
                                 "Input": target['HCOP Mouse'],
                                 "Input Type": "current symbol",
                                 "MGI Gene/Marker ID": "No mouse homolog found"}
                self.jackson_lab_data = self.jackson_lab_data.append(new_mouse_row, ignore_index=True)
        except Exception as e:
            print("Exception in JaxLab function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            self.get_jackson_lab_data(ind, target, displayname)
    
    def get_barres_lab_data(self, ind, target, displayname):
        try:
            if target['HGNC Name']:
                self.targets.at[ind, 'Barres lab human query'] = target['HGNC Name']
                # self.status.emit("Downloading data for target " + displayname + "... (Barres lab human)")
                self.emit_download_status(target, displayname, "Barres lab human")
                if target['HGNC Name'] in self.barreslab_human.index:
                    for col in self.targets.columns.tolist():
                        if col in self.barreslab_human.columns.tolist():
                            # Found a matching column; now find matching genes. 
                            # This is cunning; we look for rows in the Barres lab data with a matching gene name. 
                            pt = self.barreslab_human.ix[target['HGNC Name']][col]
                            self.targets.at[ind, col] = pt
            
            if not isnull(target['HCOP Mouse']):
                # self.targets.set_value(index=ind, col='Barres lab mouse query', value=target['MGI Symbol'])
                self.targets.at[ind, 'Barres lab mouse query'] = target['HCOP Mouse']
                # self.status.emit("Downloading data for target " + displayname + "... (Barres lab mouse)")
                self.emit_download_status(target, displayname, "Barres lab mouse")
                if target['HCOP Mouse'] in self.barreslab_mouse['Raw Data']['Gene symbol'].tolist():
                    for col in self.targets.columns.tolist():
                        if col in self.barreslab_mouse['Raw Data'].columns.tolist():
                            # Found a matching column; now find matching genes. 
                            # This is cunning; we look for rows in the Barres lab data with a matching gene name. 
                            pt = self.barreslab_mouse['Raw Data'][self.barreslab_mouse['Raw Data']['Gene symbol']
                                                                  == target['HCOP Mouse']][col].values[0]
                            self.targets.at[ind, col] = pt
        except Exception as e:
            print("Exception in Barres lab function:")
            print(e)
            print_exc()
            exit()  # There really shouldn't be any unresolvable problems here
    
    def get_pharos_data(self, ind, target, displayname):
        # Query Pharos -  an NCBI repository of drug target-related data. 
        try:
            if target['Uniprot ID'] not in [None, nan, 'nan']:
                self.emit_download_status(target, displayname, "Pharos")
                self.targets.at[ind, 'Pharos query'] = target['Uniprot ID']
                pharos_data = self.download_pharos_drug_and_ligand_count(target['Uniprot ID'])
                if pharos_data:
                    for k in pharos_data:
                        self.targets.at[ind, k] = pharos_data[k]
        except Exception as e:
            print("Exception in Pharos function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            self.get_pharos_data(ind, target, displayname)
    
    def download_pharos_drug_and_ligand_count(self, uniprot_id):
        # Deal with the complexity of digging through what a Pharos query returns.
        if uniprot_id:
            try:
                retrieval_url = "https://pharos.nih.gov/idg/api/v1/targets(" + str(uniprot_id) + ")?view=full"
                retrieval_response = get(retrieval_url)
                if retrieval_response.status_code == 404:
                    return None
                if retrieval_response.status_code != 200:
                    print("Possible Pharos problem:\nPharos retrieval query failed with code " +
                          str(retrieval_response.status_code))
                else:
                    result = loads(retrieval_response.content)
                    druggable_class = result.get('idgFamily')
                    # OK, now we're looking through a slightly deeper part of the result's structure and 
                    # counting up matching stuff... Links maybe? Yep - what I'm after is stored in a somewhat 
                    # unstructured way in [links]. 
                    ligand_count = 0
                    low_activity_ligand_count = 0
                    drug_count = 0
                    antibody_count = 0
                    literature_count = 0
                    for link in result['links']:
                        if link['kind'] == "ix.idg.models.Ligand":
                            if "Pharmalogical Action" in [i['label'] for i in link['properties']]:
                                drug_count += 1
                            if "Ligand Activity" in [i['label'] for i in link['properties']]:
                                # Implement a threshold for ligand activity here.
                                # Problem: thresholds are different for different kinds of proteins...
                                ec50 = [i['numval'] for i in link['properties'] if i['label'] in ["EC50", "IC50"]]
                                # IC50/EC50 are inverse concepts, so we should only get one or the other.
                                activity = ec50[0] if ec50 else 0
                                threshold = self.ligand_activity_thresholds.get(druggable_class) or \
                                            self.ligand_activity_thresholds['Non-IDG']
                                if activity >= threshold:
                                    ligand_count += 1
                                else:
                                    low_activity_ligand_count += 1
                    if result.get("properties"):
                        antibody_count = len([i for i in result['properties'] if
                                              i['label'] == "IDG Tools" and i['term'] == 'Antibodies'])
                    if result.get("publications"):
                        literature_count = len(result['publications'])
                    
                    returnable_data = {'Pharos ID': uniprot_id,
                                       'Druggable class': druggable_class,
                                       'Pharos': result.get('idgTDL'),
                                       'ChEMBL drug': drug_count,
                                       'ChEMBL ligand': ligand_count,
                                       'ChEMBL low-activity ligand': low_activity_ligand_count,
                                       'Antibodies': antibody_count,
                                       'Pharos literature': literature_count}
                    self.pharos_ligand_data_library[uniprot_id] = returnable_data
                    return returnable_data
            except Exception as e:
                print("Exception in Pharos drug/ligand count function:")
                print(e)
                print_exc()
                print("Retrying in 10 seconds...")
                sleep(10)
                return self.download_pharos_drug_and_ligand_count(uniprot_id)
        else:
            return None
    
    def get_disease_association_data(self, ind, target, displayname):
        # Use OpenTargets to get an absolute wealth of useful info. 
        try:
            
            if target['GeneID'] not in [None, nan, 'nan']:
                print(1)
                disease_found = False
                self.emit_download_status(target, displayname, "OpenTargets")
                self.targets.at[ind, 'OpenTargets query'] = target['GeneID']
                assocs_for_target = self.get_associations(target)
                monogenic = DataFrame(columns=self.diseasecols)
                polygenic = DataFrame(columns=self.diseasecols)
                if assocs_for_target:
                    try:
                        evidence_types = {}
                        for assoc in assocs_for_target:
                            # pprint(assoc)
                            new_dis_row = {'series': target['series'],
                                           'GeneID': target['GeneID'],
                                           'HGNC Name': target['HGNC Name'] if 'HGNC Name' in target else None,
                                           'Uniprot ID': target['Uniprot ID'] if 'Uniprot ID' in target else None,
                                           'Approved Name': target[
                                               'Approved Name'] if 'Approved Name' in target else None,
                                           'OpenTargets query': target['GeneID'], 'Source': "OpenTargets",
                                           'Disease name': assoc.get('disease').get('efo_info').get(
                                               'label') or 'Unnamed disease',
                                           'Disease ID': assoc.get('disease').get('id') or 'Unknown',
                                           'Therapeutic areas': ", ".join(
                                               assoc.get('disease').get('efo_info').get('therapeutic_area').get(
                                                   'labels')) or None,
                                           'Association score': assoc.get('association_score').get('overall') or 0,
                                           "Genetic association": assoc.get("association_score").get('datatypes').get(
                                               'genetic_association') or 0,
                                           'Direct association': True if assoc['is_direct'] else None}
                            
                            # if float(new_dis_row['Genetic association']) >= disease_association_score_threshold and \
                            #         assoc['is_direct']:
                            # self.disease_association_data = self.disease_association_data.append(new_dis_row,
                            #                                                                      ignore_index=True)
                            # This was previously part of the above conditional block
                            if match('Orphanet_', new_dis_row['Disease ID']):
                                monogenic = monogenic.append(new_dis_row, ignore_index=True)
                            else:
                                polygenic = polygenic.append(new_dis_row, ignore_index=True)
                            disease_found = True
                            
                            if evidence_types:
                                for i in assoc['evidence_count']['datatypes']:
                                    evidence_types[i] += assoc['evidence_count']['datatypes'][i]
                            else:
                                evidence_types = {i: assoc['evidence_count']['datatypes'][i] for i in
                                                  assoc['evidence_count']['datatypes']}
                            
                            # OpenTargets gives us a few stats concerning this too now!
                            if assoc.get('target'):
                                self.targets.at[ind, "OpenTargets Ensemble"] = assoc.get('target').get(
                                    'tractability').get('smallmolecule').get('ensemble')
                                self.targets.at[ind, "OpenTargets high quality compounds"] = assoc.get('target').get(
                                    'tractability').get('smallmolecule').get('high_quality_compounds')
                                self.targets.at[ind, "Small molecule genome member"] = assoc.get('target').get(
                                    'tractability').get('smallmolecule').get('small_molecule_genome_member')
                        
                        self.opentargets_datatypes[target['GeneID']] = evidence_types
                        for i in evidence_types:
                            self.targets.at[ind, sub("_", " ", i).capitalize()] = evidence_types[i]
                    except Exception as e:
                        print("Error in disease association data:")
                        print(e)
                        print_exc()
                
                # That may also be the best way to go with OMIM data. Fit it into the same table, but clearly indicate 
                # that it's from a different source. 
                # The 'evidence for target' function is even more of a goldmine. In 
                # here, you'll see Orphanet IDs, links to ChEMBL molecules, known drugs, known mechanisms of action, 
                # literature sources and other provenance info, and more! Nice. This would be good for making an 
                # 'existing drugs' table, actually. 
                print(2)
                evidence_for_target = self.get_evidence(target)
                drug2clin_evidence = []
                literature_evidence = []
                drug_found = False
                if evidence_for_target:
                    for ev in evidence_for_target:
                        if ev.get('evidence').get('drug2clinic'):
                            # This  just pulls out the bit relating to clinical trials, if there is anything.
                            try:
                                drug2clin_evidence.append(evidence_for_target['evidence']['drug2clinic'])
                            except Exception as e:
                                # This seems to often error out for no readily apparent reason. Just quietly handle it.
                                pass
                        if ev.get('drug'):
                            try:
                                # Bear in mind that many targets won't have any known drugs! 
                                new_ev_row = {'series': target['series'],
                                              'Boxed warning': False,
                                              'GeneID': target['GeneID'],
                                              'HGNC Name': target['HGNC Name'] if 'HGNC Name' in target else None,
                                              'OpenTargets query': target['GeneID'],
                                              'Disease': ev.get('disease').get('efo_info').get(
                                                  'label') or 'Unnamed disease',
                                              'Association score': ev.get('scores').get('association_score') or 0,
                                              'Drug name': ev.get('drug').get('molecule_name') or 'Unnamed molecule',
                                              'Molecule type': ev.get('drug').get('molecule_type') or 'Unknown',
                                              'Action type': ev.get('evidence').get('target2drug').get('action_type')
                                                             or 'Unknown',
                                              'Mechanism of action': ev.get('evidence').get('target2drug').
                                                                         get('mechanism_of_action')
                                                                     or 'Unknown',
                                              'Associated': ev.get('evidence').get('target2drug').get('is_associated')
                                                            or 'Unknown',
                                              'Resource score': ev.get('evidence').get('target2drug').get(
                                                  'resource_score').
                                                                    get('value') or 0,
                                              'Max clinical phase': ev.get('evidence').get('drug2clinic').
                                                                        get('clinical_trial_phase').
                                                                        get('numeric_index') or 0,
                                              'Drug ID': ev.get('drug').get('id') or 'Unnamed drug',
                                              'Withdrawn': None,
                                              'Withdrawn reason': ev.get('drug').get('withdrawn_reason'),
                                              'Withdrawn country': ev.get('drug').get('withdrawn_country'),
                                              'Withdrawn year': ev.get('drug').get('withdrawn_year')
                                              }
                                if new_ev_row['Withdrawn reason']:
                                    new_ev_row['Withdrawn'] = True
                                # Has this drug been withdrawn? Check the WITHDRAWN database.
                                if new_ev_row.get('Drug name'):
                                    if str(new_ev_row['Drug name']).capitalize() in self.withdrawn_withdrawn[
                                        'DRUG_NAME'].tolist():
                                        new_ev_row['Withdrawn'] = True
                                # To prevent errors, we'll have to do some catches for missing values here.
                                if float(new_ev_row['Association score']) >= disease_association_score_threshold:
                                    self.existing_drug_data = self.existing_drug_data.append(new_ev_row,
                                                                                             ignore_index=True)
                                    drug_found = True
                                    # Set these too - they appear in the SM Druggability sheet.
                                    # It's the same info, but in a potentially more convenient arrangement.
                                    if new_ev_row['Withdrawn']:
                                        self.targets.at[ind, 'Has withdrawn drug'] = True
                                        druglist = []
                                        if target['Withdrawn drug list']:
                                            druglist = target['Withdrawn drug list'].split(", ")
                                        druglist.append(new_ev_row['Drug name'])
                                        self.targets.at[ind, 'Withdrawn drug list'] = ", ".join(druglist)
                            except Exception as e:
                                print("Error in disease association data:")
                                print(e)
                                print_exc()
                        # This is also where we can get Orphanet IDs, so we should think about adding them to the mix 
                        # somehow. Note that a lot of entries in this list won't have Orphanet IDs (because they don't 
                        # exist), but some will. Orphanet IDs should be added to disease_association_data.
                        # Source = OrphaNet, Disease ID = orphanet ID if 'name' not in ev['disease']: pprint(ev) input() 
                        # if match("Orphanet_", ev['disease']['id']):
                        #     try:
                        #         new_dis_row = {'series': target['series'], 'GeneID': target['GeneID'],
                        #                        'HGNC Name': target['HGNC Name'],
                        #                        'Approved Name': target['Approved Name'],
                        #                        'OpenTargets query': target['GeneID'], 'Source': "Orphanet",
                        #                        'Disease name': ev.get('disease').get('efo_info').get(
                        #                            'label') or 'Unnamed disease',
                        #                        'Disease ID': ev.get('disease').get('id') or 'Unnamed disease',
                        #                        'Association score': ev.get('scores').get('association_score') or 0}
                        #         
                        #         if float(new_dis_row['Association score']) >= disease_association_score_threshold:
                        #             # self.disease_association_data = self.disease_association_data.append(new_dis_row,
                        #             #                                                                      ignore_index=True)
                        #             if match('Orphanet_', new_dis_row['Disease ID']):
                        #                 monogenic = monogenic.append(new_dis_row, ignore_index=True)
                        #             else:
                        #                 polygenic = polygenic.append(new_dis_row, ignore_index=True)
                        #             disease_found = True
                        #     except Exception as e:
                        #         print("Error in disease association data:")
                        #         print(e)
                        #         print_exc()
                        # There's yet more we can get from here. We can get papers too.
                        if ev.get("literature"):
                            # pprint(js['literature'])
                            if ev['literature'].get('title'):
                                authors = ['Unknown']
                                if ev['literature'].get('authors'):
                                    authors = [i.get('short_name') for i in ev['literature']['authors'] if
                                               'short_name' in i]
                                ln = {
                                    'series': target['series'],
                                    'HGNC Name': target['HGNC Name'] if 'HGNC Name' in target else None,
                                    'GeneID': target['GeneID'] if 'GeneID' in target else None,
                                    'Uniprot ID': target['Uniprot ID'] if 'Uniprot ID' in target else None,
                                    "Title": ev['literature']['title'],
                                    "Authors": ", ".join(authors),
                                    "Date": ev['literature'].get("date") or "",
                                    "Journal ref": ev['literature'].get("journal_reference") or "",
                                    "DOI": ev['literature'].get("doi") or ""
                                }
                                literature_evidence.append(ln)
                
                if drug2clin_evidence:
                    self.clinical_trial_data[target['GeneID']] = drug2clin_evidence
                    for dt in drug2clin_evidence:  # dt = data type
                        colname = sub("_", " ", dt).capitalize()
                        self.targets.at[ind, colname] = drug2clin_evidence[dt]
                
                if literature_evidence and self.get_literature:
                    self.literature = concat([self.literature, DataFrame(literature_evidence)], sort=False)
                
                # Before actually committing disease associations to more permanent storage, we need a bit of clean-up 
                # on them. Sort by association score and remove duplicates.
                monogenic = monogenic.sort_values(by=['Genetic association'])
                monogenic = monogenic.drop_duplicates(['HGNC Name', 'GeneID', 'Uniprot ID', 'Disease ID'])
                for i, row in monogenic.iterrows():
                    self.monogenic_disease_association_data = self.monogenic_disease_association_data.append(
                        row, ignore_index=True)
                
                polygenic = polygenic.sort_values(by=['Genetic association'])
                polygenic = polygenic.drop_duplicates(['HGNC Name', 'GeneID', 'Uniprot ID', 'Disease ID'])
                for i, row in polygenic.iterrows():
                    self.polygenic_disease_association_data = self.polygenic_disease_association_data.append(
                        row, ignore_index=True)
                
                if not drug_found:
                    # If no drugs found, leave a row in self.existing_drug_data explicitly stating it.
                    new_drug_row = {
                        'series': target['series'],
                        'GeneID': target['GeneID'],
                        'HGNC Name': target['HGNC Name'],
                        'OpenTargets query': target['GeneID'],
                        'Drug name': "No existing drugs found",
                        'Boxed warning': False
                    }
                    self.existing_drug_data = self.existing_drug_data.append(new_drug_row, ignore_index=True)
                
                if not disease_found:
                    # If no diseases found, leave a row explicitly stating it.
                    new_dis_row = {
                        'series': target['series'],
                        'GeneID': target['GeneID'],
                        'HGNC Name': target['HGNC Name'],
                        'Approved Name': target['Approved Name'],
                        'OpenTargets query': target['GeneID'],
                        'Source': "No disease associations found"
                    }
                    self.monogenic_disease_association_data = self.monogenic_disease_association_data.append(
                        new_dis_row, ignore_index=True)
                    self.polygenic_disease_association_data = self.polygenic_disease_association_data.append(
                        new_dis_row, ignore_index=True)
        except Exception as e:
            print("Exception in OpenTargets function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            self.get_disease_association_data(ind, target, displayname)
    
    def get_associations(self, target):
        # Make one type of OpenTargets query. 
        # (This function exists for error-handling purposes).
        try:
            assocs = self.ot.get_associations_for_target(target['GeneID'])
        except Exception as e:
            print("Exception in Disease Association function, pt1:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            assocs = self.get_associations(target)
        return assocs
    
    def get_evidence(self, target):
        # Make the other type of OpenTargets query. 
        # (This function exists for error-handling purposes).
        try:
            evidence = self.ot.get_evidence_for_target(target['GeneID'])
        except Exception as e:
            print("Exception in Disease Association function, pt2:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            evidence = self.get_evidence(target)
        return evidence
    
    def get_drug_safety_data(self, ind, target, displayname):
        # Source of drug names: self.existing_drug_data (see self.get_disease_association_data())
        # I did it this way because we often get a LOT of duplicate drug names from OpenTargets. This is much more 
        # efficient.
        try:
            if target['HGNC Name'] in self.existing_drug_data['HGNC Name'].tolist():
                # self.status.emit("Downloading data for target " + displayname + "... (Drug safety warnings)")
                self.emit_download_status(target, displayname, "FDA safety warnings")
                drugs = self.existing_drug_data[self.existing_drug_data['HGNC Name'] == target['HGNC Name']][
                    'Drug name'].tolist()
                # This list is likely to have duplicate entries.
                drugs = [i for i in uniq(drugs) if i != 'Unnamed molecule']
                print("got some drugs for this target")
                drugs_with_warnings = []
                for d in drugs:
                    print("    " + d)
                    if self.api_keys['FDA_key'] == "FDAkey":
                        url = "https://api.fda.gov/drug/label.json?search=substance_name:" + d
                    else:
                        url = "https://api.fda.gov/drug/label.json?api_key=" + self.api_keys[
                            'FDA_key'] + "&search=substance_name:" + d
                    response = get(url)
                    if response.status_code == 200:
                        data = loads(response.text)
                        if data.get("error"):
                            pprint(data)  # Probably can comment this out later, it'll be mostly NOT FOUND, which is OK
                        elif data.get("results"):
                            # This is where we deal with actual results.
                            if data['results'][0].get('boxed_warning'):
                                drugs_with_warnings.append(d)
                                print("        HAS A WARNING")
                        else:
                            pprint(data)  # This is likely to be where we see real trouble.
                    elif response.status_code == 404:
                        pass
                    else:
                        print("It's all gone wrong in the FDA query: code " + str(response.status_code))
                # Good - now modify the existing dataframe to hold that info.
                if drugs_with_warnings:
                    self.existing_drug_data['Boxed warning'][
                        self.existing_drug_data['Drug name'].isin(drugs_with_warnings)] = True
        except Exception as e:
            print("Exception in FDA drug safety function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            self.get_drug_safety_data(ind, target, displayname)
    
    def get_gnomad_data(self, ind, target, displayname):
        # Collect a couple of stats from the gnomAD database.
        try:
            if target['HGNC Name']:
                if target['HGNC Name'] in self.gnomad['gene'].tolist():
                    # self.status.emit("Downloading data for target " + displayname + "... (Drug safety warnings)")
                    self.emit_download_status(target, displayname, "gnomAD")
                    gn = self.gnomad[self.gnomad['gene'] == target['HGNC Name']].sort_values(by=['pLI'],
                                                                                             ascending=False).iloc[0]
                    self.targets.at[ind, 'gnomAD pLI'] = gn['pLI']
                    self.targets.at[ind, 'gnomAD o/e'] = gn['oe_lof']
        except Exception as e:
            print("Exception in gnomAD data function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            self.get_drug_safety_data(ind, target, displayname)
    
    def get_gwas_data(self, ind, target, displayname):
        # Look up gene in two different GWAS experiment libraries.
        try:
            if target['HGNC Name'] not in [None, nan, 'nan']:
                # self.status.emit("Downloading data for target " + displayname + "... (GWAS)")
                self.emit_download_status(target, displayname, "GWAS")
                self.targets.at[ind, 'GWAS query'] = target['HGNC Name']
                gwas_count = 0
                # url = "http://www.gwascentral.org/studies?q=" + target['HGNC Name'] + "&t=6&format=json"
                # print(url)
                # response = get(url)
                # if response.status_code != 200:
                #     pass
                # else:
                #     for gwas in loads(response.text):
                #         # Bear in mind that many - possibly most - targets won't have any known GWAS findings! 
                #         new_gwas_row = {
                #             'series': target['series'],
                #             'GeneID': target['GeneID'],
                #             'HGNC Name': target['HGNC Name'],
                #             'Source': 'GWAS Central',
                #             'Study name': gwas['name'],
                #             'ID': gwas['identifier'],
                #             'Phenotypes': ", ".join(gwas['phenotypes']),
                #             'Highest P-value': gwas['highest_pvalue'],
                #             'Num. markers': gwas['number_of_markers'],
                #             'link': gwas['link']
                #         }
                #         self.gwas_data = self.gwas_data.append(new_gwas_row, ignore_index=True)
                #         gwas_count += 1
                # Disabled because GWAS Central have stopped allowing bulk querying via API...
                # Since we're doing this anyway, we may as well search the GWAS Catalog data that I opened earlier.
                for x, gwas in self.gwas_catalog[
                            self.gwas_catalog['MAPPED_GENE'] == target['HGNC Name'].upper()].iterrows():
                    new_gwas_row = {
                        'series': target['series'],
                        'GeneID': target['GeneID'],
                        'HGNC Name': target['HGNC Name'],
                        'Source': 'GWAS Catalog',
                        'Study name': gwas['STUDY'],
                        'ID': gwas['STUDY ACCESSION'],
                        'Phenotypes': gwas['MAPPED_TRAIT'],
                        'Highest P-value': gwas['P-VALUE'],
                        'link': gwas['LINK']
                    }
                    self.gwas_data = self.gwas_data.append(new_gwas_row, ignore_index=True)
                    gwas_count += 1
                
                if gwas_count == 0:
                    new_gwas_row = {
                        'series': target['series'],
                        'GeneID': target['GeneID'],
                        'HGNC Name': target['HGNC Name'],
                        'Source': 'No GWAS publications found'
                    }
                    self.gwas_data = self.gwas_data.append(new_gwas_row, ignore_index=True)
        except Exception as e:
            print("Exception in GWAS function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            self.get_gwas_data(ind, target, displayname)
    
    def get_antibodyability(self, ind, target, displayname):
        # AB-ability is a hodgepodge of various bits of data, many of which we've actually already got (mostly from
        # the HPA). Here we'll get the remaining stuff from several sources. 
        try:
            self.emit_download_status(target, displayname, "AB-ability")
            # Is this target a member of the surfaceome? (Human or mouse - either is interesting)
            # (The space in "Protein " is intentional; that's what's in the header)
            # If so, what level of confidence are we looking at?
            if target['Uniprot ID'] in self.surfaceome_inclusion["Table A (human)"]['Protein '].values:
                self.targets.at[ind, 'Surfaceome membership (human)'] = True
                # Take subset of data matching uniprot ID, get first confidence value in list (they're all the same)
                confidence = self.surfaceome_inclusion["Table A (human)"][
                    self.surfaceome_inclusion["Table A (human)"]['Protein '] == target['Uniprot ID']][
                    'CSPA category'].values.tolist()[0]
                self.targets.at[ind, 'Surfaceome membership confidence (human)'] = confidence
            if target['Mouse Uniprot ID']:
                if any(x in self.surfaceome_inclusion["Table B (mouse)"]['Protein '].values
                       for x in target['Mouse Uniprot ID']):
                    self.targets.at[ind, 'Surfaceome membership (mouse)'] = True
                    confidence = self.surfaceome_inclusion["Table B (mouse)"][
                        self.surfaceome_inclusion["Table B (mouse)"]['Protein '].isin(target['Mouse Uniprot ID'])][
                        'CSPA category'].values.tolist()[0]
                    self.targets.at[ind, 'Surfaceome membership confidence (mouse)'] = confidence
        except Exception as e:
            print("Exception in AB-ability function:")
            print(e)
            print_exc()
    
    def get_pdb_structures(self, ind, target, displayname):
        # Query the Protein Data Bank for records of 3D structure models.
        try:
            if target['Uniprot ID']:
                self.emit_download_status(target, displayname, "PDB structures")
                pdb_ids = self.pdb_to_uniprot[self.pdb_to_uniprot["UniProt ID"] == target['Uniprot ID']][
                    'PDB ID'].tolist()
                described_pdbs = []
                for i in pdb_ids:
                    # Occasionally this fails, apparently due to a bug in pypdb. This is beyond my control, so we'll 
                    # have to try-except our way around it.
                    try:
                        described_pdbs.append(describe_pdb(i))
                    except Exception:
                        described_pdbs.append({'structureId': i})
                df = DataFrame(described_pdbs)
                
                # Now add back in the usual IDs, and we're almost there.
                idcols = ['series', 'GeneID', 'HGNC Name', 'Uniprot ID']
                for col in idcols:
                    if target.get(col):
                        df[col] = target[col]
                # Swap those IDs to be at the front
                df = df[df.columns.tolist()[-4:] + df.columns.tolist()[:-4]]
                self.pdb_data = concat([self.pdb_data, df], ignore_index=True, sort=False)
                # Got to force these columns to stay in order
                self.pdb_data = self.pdb_data[idcols + [i for i in self.pdb_data.columns.tolist() if i not in idcols]]
        except Exception as e:
            print("Exception in PDB structure function:")
            print(e)
            print_exc()
            exit()
    
    def get_antitargets(self, ind, target, displayname):
        if target['Uniprot ID']:
            self.emit_download_status(target, displayname, "Antitargets")
            i = target['Uniprot ID']
            homologs = self.get_homologs(i)
            if homologs:
                for j in [str(k) for k in homologs]:
                    if i != j:
                        identity = (self.identity_scores[i][j] + self.identity_scores[j][i]) / 2
                        d = target[self.antitargetcols[:4]]
                        d['Antitarget Uniprot ID'] = j
                        d['Antitarget gene name'] = self.get_gene_name_for_uniprot_id(j)
                        d['Protein seq. identity (%)'] = identity
                        if identity >= self.api_keys['seq_similarity_threshold']:
                            self.antitargets_data = self.antitargets_data.append(d, ignore_index=True)
            else:
                d = target[self.antitargetcols[:4]]
                d['Antitarget Uniprot ID'] = None
                d['Antitarget gene name'] = None
                d['Protein seq. identity (%)'] = 0
                self.antitargets_data = self.antitargets_data.append(d, ignore_index=True)
    
    def get_protein_protein_interactions(self, ind, target, displayname):
        # Query BioGRID for protein-protein interactions.
        try:
            if target['HGNC Name'] and self.get_ppi:
                self.emit_download_status(target, displayname, "protein-protein interactions")
                url = "https://webservice.thebiogrid.org/interactions/?accesskey=" + self.api_keys[
                    'biogrid_key'] + "&taxId=9606&includeHeader=true&interSpeciesExcluded=true&geneList=" + target[
                          'HGNC Name']
                df = read_table(url)  # Pandas is so great
                # That gives us a lot of columns, but we don't need all of them.
                dropcols = ['Entrez Gene Interactor A',
                            'Entrez Gene Interactor B',
                            'BioGRID ID Interactor A',
                            'BioGRID ID Interactor B',
                            'Systematic Name Interactor A',
                            'Systematic Name Interactor B', ]
                for i in dropcols:
                    df = df.drop(i, axis=1)
                # Now add back in the usual IDs, and we're almost there.
                idcols = ['series', 'GeneID', 'HGNC Name', 'Uniprot ID']
                for col in idcols:
                    if target.get(col):
                        df[col] = target[col]
                    else:
                        df[col] = None
                # Swap those IDs to be at the front
                df = df[df.columns.tolist()[-4:] + df.columns.tolist()[:-4]]
                self.ppi_data = concat([self.ppi_data, df], ignore_index=True, sort=False)
                # Got to force these columns to stay in order
                self.ppi_data = self.ppi_data[idcols + [i for i in self.ppi_data.columns.tolist() if i not in idcols]]
        except Exception as e:
            print("Exception in protein-protein interaction function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            self.get_protein_protein_interactions(ind, target, displayname)
    
    def get_drug_toxicology(self, ind, target, displayname):
        # Query T3DB toxicology data.
        # This is some great stuff and I'm really pleased to have found it. 
        try:
            if target['Uniprot ID']:
                df = self.t3db_drug_moas[self.t3db_drug_moas['Target UniProt ID'] == target['Uniprot ID']].copy()
                # Easiest thing in the world. Now add back in the usual IDs, and we're almost there.
                idcols = ['series', 'GeneID', 'HGNC Name']
                for col in idcols:
                    if target.get(col):
                        df[col] = target[col]
                    else:
                        df[col] = None
                # Swap those IDs to be at the front
                df = df[df.columns.tolist()[-3:] + df.columns.tolist()[:-3]]
                self.drug_tox_data = concat([self.drug_tox_data, df], ignore_index=True, sort=False)
                # Got to force these columns to stay in order
                self.drug_tox_data = self.drug_tox_data[idcols + [i for i in df.columns.tolist() if i not in idcols]]
        except Exception as e:
            print("Exception in PDB structure function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            self.get_pdb_structures(ind, target, displayname)
    
    def gather_output_sheets(self):
        # Put some sensible output together.
        # Instead of writing that all to a single sheet, try picking out relevant parts and writing those to 
        # different sheets. Do it at the end here to minimise potential for error and confusion. 
        # Input data sheet
        # All that really matters here is that these few columns go first, in this order:
        sheets_by_name = OrderedDict()
        core_cols = ['series', 'HGNC Name', 'GeneID', 'Uniprot ID']
        self.input_data = self.input_data[[i for i in core_cols if i in self.input_data.columns.tolist()] +
                                          [i for i in self.input_data.columns.tolist() if i not in core_cols]]
        sheets_by_name['External data input'] = self.input_data
        # Basic information sheet 
        print("Prepping output sheets")
        print("basic")
        sheets_by_name['Basic information'] = self.targets[self.basicinfocols].copy()
        # Buckets sheet
        print("buckets")
        sheets_by_name['Buckets'] = self.targets[self.bucketcols].copy()
        sheets_by_name['Bucket guide'] = self.bucketing_guide.copy()
        # Safety risk sheet
        print("safety")
        sheets_by_name['GTEX'] = self.targets[self.gtexcols].copy()
        print("hpa expression")
        sheets_by_name['HPA expression'] = self.hpa_expression_data.copy()
        # Barres mouse sheet
        print("barres mouse")
        sheets_by_name['Barres mouse'] = self.targets[self.barresmousecols]
        # Barres human sheet
        print("barres human")
        sheets_by_name['Barres human'] = self.targets[self.barreshumancols]
        # HPA enrichment sheet
        print("hpa enrichment")
        sheets_by_name['HPA Enrichment'] = self.targets[
            self.hpaenrichmentcols + ["HPA " + i for i in hpa_all_tissues]]
        # Risk factors sheet
        print("risk")
        sheets_by_name['Risk factors'] = self.targets[self.riskfactorcols].copy()
        # WikiPathways sheet
        print("pathways")
        sheets_by_name['Pathways'] = self.pathways.copy()
        # Disease associations sheet
        # Actually, this is WAY simpler:
        # Columns specified in order to set their printed order
        print("disease")
        if self.monogenic_disease_association_data.columns.tolist():
            sheets_by_name['Rare disease associations'] = self.monogenic_disease_association_data[self.diseasecols]
        if self.polygenic_disease_association_data.columns.tolist():
            sheets_by_name['Common disease associations'] = self.polygenic_disease_association_data[self.diseasecols]
        # MGI mouse sheet
        print("jaxlab")
        sheets_by_name['MGI mouse'] = self.jackson_lab_data[self.jacksonlabcols]
        # CanSAR sheet
        print("cansar")
        sheets_by_name['CanSAR'] = self.targets[self.cansarcols]
        # PDB data sheet
        print("pdb")
        sheets_by_name['Protein structures'] = self.pdb_data.copy()  # No mods needed - it's grown in a safe manner.
        # Antitargets sheet
        print("antitargets")
        sheets_by_name['Antitargets'] = self.antitargets_data[self.antitargetcols]
        # Core fitness gene sheet
        print("core fitness")
        sheets_by_name['Core fitness'] = self.targets[self.corefitnesscols]
        # Druggability sheet
        print("sm druggability")
        sheets_by_name['SM Druggability'] = self.targets[self.druggabilitycols].copy()
        # AB-ability sheet
        print("antibodyability")
        sheets_by_name['AB-ability'] = self.targets[self.antibodyabilitycols].copy()
        # Feasibility sheet
        print("feasibility")
        sheets_by_name['Feasibility'] = self.targets[self.feasibilitycols].copy()
        # Existing drugs sheet
        # Columns specified in order to set their printed order
        print("drugs")
        if self.existing_drug_data.columns.tolist():
            # This sheet appears to have duplicate rows sometimes in the final rendition. 
            sheets_by_name['Existing drugs'] = self.existing_drug_data[self.existingdrugscols].drop_duplicates()
        # Drug toxicology sheet
        sheets_by_name['Drug toxicology'] = self.drug_tox_data.copy()  # No mods - grown in a safe manner.
        # Pharos sheet
        print("pharos")
        sheets_by_name['Pharos'] = self.targets[self.pharoscols].copy()
        # GWAS sheet
        # Columns specified in order to set their printed order
        print("gwas")
        if self.gwas_data.columns.tolist():
            sheets_by_name['GWAS'] = self.gwas_data[self.gwascols]
        # Protein-protein interaction data sheet
        # This gives us a pretty huge output sheet. I may switch this off in the future. 
        if self.get_ppi:
            print("protein-protein interactions")
            sheets_by_name['Protein-protein interactions'] = self.ppi_data.copy()  # No mods - grown in a safe manner.
        # literature
        if self.get_literature:
            print("literature")
            sheets_by_name['Literature'] = self.literature[self.literaturecols].copy()
        
        # Covert numbers to numeric?
        # for i in sheets_by_name:
        #     sheets_by_name[i] = sheets_by_name[i].apply(to_numeric, errors='ignore')
        
        return sheets_by_name
    
    def get_ensembl_id_for_gene_name(self, gene_name):
        # Use the lookup table first; should that fail, try querying genenames.
        try:
            df = self.id_lookup_table[self.id_lookup_table['Approved Symbol'] == gene_name]
            if len(df):
                ensembl_id = str(df['Ensembl Gene ID'].tolist()[0]).strip()
                # There may be cases where a record exists but no actual value is supplied, or the value is incorrect in 
                # some way. Check for it.
                if match("ENSG[0-9]{11}", ensembl_id):
                    return ensembl_id
        except Exception as e:
            print("Exception in HGNC Name -> Ensembl ID function:")
            print(e)
            print_exc()
            exit()
        # Should that fail, go and look it up from Ensembl's servers
        try:
            url = "http://rest.genenames.org/fetch/symbol/" + gene_name
            print("Looking up ensembl ID for gene name " + gene_name)
            print("url = " + url)
            response = get(url)
            if response.status_code != 200:
                # self.warnings.emit(
                #     "GeneNames query gene name -> Ensembl ID query failed with code " + str(response.status_code))
                if response.status_code == 500:
                    sleep(10)
                    ensembl_id = self.get_ensembl_id_for_gene_name(gene_name)
                    return ensembl_id
            else:
                data = xmltodict.parse(response.content)
                # Does the following bit of xml exist?
                if int(data.get('response').get('result').get('@numFound')) > 0:
                    ensembl_bits = [i for i in data['response']['result']['doc']['str'] if
                                    i['@name'] == "ensembl_gene_id"]
                    # This can, occasionally, be empty...
                    if len(ensembl_bits) > 0:
                        ensembl_id = ensembl_bits[0]['#text']
                        if ensembl_id is nan:
                            ensembl_id = None
                        return ensembl_id
                    else:
                        return None
                        # return ensembl_id, None
                else:
                    # In this case, no matching record in genenames has been found. 
                    # This could be because the user has supplied a synonymous or outdated gene symbol. Try to use 
                    # BioMart to resolve that.
                    # No dice. I'll come back to this later, but for now, the program should refuse to work here until a
                    # correct name is supplied. 
                    # exit("ERROR: Unable to find records for gene symbol " + gene_name + "!")
                    return None
        except Exception as e:
            print("Exception in HGNC Name -> Ensembl ID function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            return self.get_ensembl_id_for_gene_name(gene_name)
    
    def get_uniprot_id_for_gene_name(self, gene_name):
        # Use the lookup table first; should that fail, try querying genenames.
        try:
            df = self.id_lookup_table[self.id_lookup_table['Approved Symbol'] == gene_name]
            if len(df):
                uniprot_id = str(df['UniProt ID(supplied by UniProt)'].tolist()[0]).strip()
                # There may be cases where a record exists but no actual value is supplied, or the value is incorrect in 
                # some way. Check for it.
                if match("[A-Z0-9]+", uniprot_id):
                    return uniprot_id
        except Exception as e:
            print("Exception in HGNC Name -> UniProt ID function:")
            print(e)
            print_exc()
            exit()
        # Should that fail, go and look it up from UniProt's servers
        try:
            url = "http://www.uniprot.org/uniprot/?query=gene:" + gene_name + "+AND+organism:" + \
                  self.organism + "&sort=score&format=tab"
            print("Looking up uniprot ID for gene name " + gene_name)
            print("url = " + url)
            response = get(url)
            uniprot_data = response.text
            # if response.status_code != 200:
            #     self.warnings.emit(
            #         "UniProt query gene name -> UniProt ID failed with code " + str(response.status_code))
            # else:
            if response.status_code in [200]:
                if isinstance(uniprot_data, str):
                    if len(uniprot_data) > 0:
                        uniprot = read_table(StringIO(uniprot_data))
                        reviewed = uniprot[uniprot['Status'] == 'reviewed']
                        if len(reviewed) > 0:
                            # We want better handling of cases where there is more than one available Uniprot ID.
                            uniprot_id = str(reviewed['Entry'].tolist()[0])
                            return uniprot_id
                else:
                    print(uniprot_data)
                return None
        except Exception as e:
            print("Exception in HGNC name -> UniProt ID function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            return self.get_uniprot_id_for_gene_name(gene_name)
    
    def get_gene_name_for_ensembl_id(self, ensembl_id):
        # Use the lookup table first; should that fail, try querying genenames.
        try:
            df = self.id_lookup_table[self.id_lookup_table['Ensembl Gene ID'] == ensembl_id]
            if len(df):
                gene_name = str(df['Approved Symbol'].tolist()[0]).strip()
                # There may be cases where a record exists but no actual value is supplied, or the value is incorrect in 
                # some way. Check for it.
                if match("[a-zA-Z0-9]+", gene_name):
                    return gene_name
        except Exception as e:
            print("Exception in Ensembl ID -> HGNC name function:")
            print(e)
            print_exc()
            exit()
        # Should that fail, go and look it up from Ensembl's servers
        try:
            url = "http://rest.ensembl.org/xrefs/id/" + ensembl_id + "?content-type=application/json"
            print("Looking up gene name for ensembl ID " + ensembl_id)
            print("url = " + url)
            response = get(url)
            # if response.status_code not in [200, 400]:
            #     self.warnings.emit(
            #         "Ensembl query Ensembl ID -> gene name failed with code " + str(response.status_code))
            # else:
            if response.status_code in [200, 400]:
                ensembl_json = loads(response.text)
                print(ensembl_json)
                if 'error' in ensembl_json:
                    return None
                else:
                    for i in ensembl_json:
                        if i['db_display_name'] == "NCBI gene":
                            return i['display_id']
                    return None
        except Exception as e:
            print("Exception in HGNC name -> Ensembl ID function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            return self.get_gene_name_for_ensembl_id(ensembl_id)
    
    def get_uniprot_id_for_ensembl_id(self, ensembl_id):
        # Use the lookup table first; should that fail, try querying genenames.
        try:
            df = self.id_lookup_table[self.id_lookup_table['Ensembl Gene ID'] == ensembl_id]
            if len(df):
                uniprot_id = str(df['UniProt ID(supplied by UniProt)'].tolist()[0]).strip()
                # There may be cases where a record exists but no actual value is supplied, or the value is incorrect in 
                # some way. Check for it.
                if match("[A-Z0-9]+", uniprot_id):
                    return uniprot_id
        except Exception as e:
            print("Exception in Ensembl ID -> UniProt ID function:")
            print(e)
            print_exc()
            exit()
        # Should that fail, go and look it up from UniProt's servers
        try:
            url = "http://rest.ensembl.org/xrefs/id/" + ensembl_id + "?content-type=application/json"
            print("Looking up uniprot ID for ensembl ID " + ensembl_id)
            print("url = " + url)
            response = get(url)
            # if response.status_code not in [200, 400]:
            #     self.warnings.emit(
            #         "Ensembl query Ensembl ID -> UniProt ID failed with code " + str(response.status_code))
            # else:
            if response.status_code in [200, 400]:
                ensembl_data = response.text
                ensembl_json = loads(ensembl_data)
                if 'error' in ensembl_json:
                    return None
                else:
                    for i in ensembl_json:
                        if i['db_display_name'] == "UniProtKB Gene Name":
                            return i['primary_id']
                    return None
        except Exception as e:
            print("Exception in Ensembl ID -> UniProt ID function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            return self.get_uniprot_id_for_ensembl_id(ensembl_id)
    
    def get_ensembl_id_for_uniprot_id(self, uniprot_id):
        # Bit tricky, this one. 
        # I didn't find a way to do it directly by query, so if we can't find it in the lookup table, we can't have it.
        # As such, I've avoided it where possible.
        try:
            df = self.id_lookup_table[self.id_lookup_table['UniProt ID(supplied by UniProt)'] == uniprot_id]
            if len(df):
                ensembl_id = str(df['Ensembl Gene ID'].tolist()[0]).strip()
                # There may be cases where a record exists but no actual value is supplied, or the value is incorrect in 
                # some way. Check for it.
                if match("ENSG[0-9]{11}", ensembl_id):
                    return ensembl_id
        except Exception as e:
            print("Exception in UniProt ID -> Ensembl ID function:")
            print(e)
            print_exc()
            exit()
    
    def get_gene_name_for_uniprot_id(self, uniprot_id):
        # Use the lookup table first; should that fail, try querying genenames.
        try:
            df = self.id_lookup_table[self.id_lookup_table['UniProt ID(supplied by UniProt)'] == uniprot_id]
            if len(df):
                gene_name = str(df['Approved Symbol'].tolist()[0]).strip()
                # There may be cases where a record exists but no actual value is supplied, or the value is incorrect in 
                # some way. Check for it.
                if match("[a-zA-Z0-9]+", gene_name):
                    return gene_name
        except Exception as e:
            print("Exception in HGNC Name -> UniProt ID function:")
            print(e)
            print_exc()
            exit()
        # Should that fail, go and look it up from UniProt's servers
        try:
            url = "http://www.uniprot.org/uniprot/?query=accession:" + uniprot_id + "+AND+organism:" + \
                  self.organism + "&sort=score&format=tab"
            print("Looking up gene name for uniprot ID " + uniprot_id)
            print("url = " + url)
            response = get(url)
            # if response.status_code != 200:
            #     self.warnings.emit(
            #         "UniProt query UniProt ID -> gene name failed with code " + str(response.status_code))
            # else:
            if response.status_code == 200:
                uniprot_data = response.text
                if uniprot_data:
                    uniprot = read_table(StringIO(uniprot_data))
                    names = []
                    if len(uniprot['Entry name'].dropna()):
                        names = [sub("_HUMAN", '', i) for i in uniprot['Entry name'].tolist()[0].split(' ')]
                    elif len(uniprot['Gene names'].dropna()):
                        names = uniprot['Gene names'].tolist()[0].split(' ')
                    if names:
                        return names[0]
                    else:
                        return uniprot_id  # Better than nothing...
        except Exception as e:
            print("Exception in HGNC Name -> UniProt ID function:")
            print(e)
            print_exc()
            print("Retrying in 10 seconds...")
            sleep(10)
            return self.get_gene_name_for_uniprot_id(uniprot_id)
    
    def set_target_list(self, targets):
        # Note the copy - this should prevent any funny business
        self.input_data = targets.copy()
    
    def bucketing(self):
        # When all downloaded information is available, we can look at assigning targets to buckets.
        print("BUCKETING")
        # Assign buckets, arrange it into a new sheet, and return a full new excel file-writable object.
        try:
            for ind, target in self.targets.iterrows():
                print(ind)
                self.emit_bucket_status(target)
                safety = self.safety_bucket(target)
                self.targets.at[ind, 'Safety bucket'] = safety
                druggability = self.druggability_bucket(target)
                self.targets.at[ind, 'SM Druggability bucket'] = druggability
                feasibility = self.feasibility_bucket(target)
                self.targets.at[ind, 'Feasibility bucket'] = feasibility
                antibodyability = self.antibodyability_bucket(target)
                self.targets.at[ind, 'AB-ability bucket'] = antibodyability
                new_modality = self.new_modality_bucket(target, druggability, antibodyability)
                self.targets.at[ind, 'New modality bucket'] = new_modality
                tissue_engagement = self.tissue_engagement_bucket(target)
                self.targets.at[ind, 'Tissue engagement bucket'] = tissue_engagement
        except Exception as e:
            print("Exception in bucketing function:")
            print(e)
            print_exc()
            exit()
    
    def safety_bucket(self, target):
        # Unlike the druggability bucket below, this is better visualised as a waterfall than a decision tree. 
        # Well, not quite. I think it can be visualised as picking the lowest result from a set of sub-decisions.
        # Though if we reach a level 2, we can potentially push it up to a level 1 if additional considerations are met.
        
        def has_insufficient_info():
            if not target['HGNC Name'] or not target['Uniprot ID']:
                return True
                # jaxlab = self.jackson_lab_data[self.jackson_lab_data['HGNC Name'] == target['HGNC Name']]
                # if "No phenotypes found" in jaxlab['MGI Gene/Marker ID'].tolist():
                #     return True
        
        def boxed_warnings():
            drugs = self.existing_drug_data[self.existing_drug_data['HGNC Name'] == target['HGNC Name']]
            drugs_with_warnings = drugs[drugs['Boxed warning'] == True]
            if len(drugs_with_warnings.index):
                return True
        
        def on_target_withdrawn_drug():
            # Check if the target has known drugs that have on-target effects (interactions) and has boxed warnings.
            # Upper-case all of their names, because some sources are inconsistently formatted.
            drugs = self.existing_drug_data[self.existing_drug_data['HGNC Name'] == target['HGNC Name']]
            withdrawn_drugs_for_this_target = drugs[drugs['Withdrawn'].notnull()]
            # If we have anything directly marked as a withdrawn drug, return a result immediately.
            if len(withdrawn_drugs_for_this_target.index):
                return True
            
            # Let's also get a list of drugs which are on record as having toxicity issues (t3db). 
            toxx = self.drug_tox_data[self.drug_tox_data['HGNC Name'] == target['HGNC Name']]
            if len(toxx.index):
                return True
        
        def on_target_boxed_warnings():
            # Check if the target has known drugs that have on-target effects (interactions) and has boxed warnings.
            # Upper-case all of their names, because some sources are inconsistently formatted.
            drugs = self.existing_drug_data[self.existing_drug_data['HGNC Name'] == target['HGNC Name']]
            drugs_with_warnings = drugs[drugs['Boxed warning'] == True]
            drugs_with_warnings_names = [i.upper() for i in drugs_with_warnings['Drug name'].tolist()]
            # Look up names of interacting drugs from DGIdb
            interactions = self.dgidb_interactions[self.dgidb_interactions['gene_name'] == target['HGNC Name']]
            interactions_names = [str(i).upper() for i in interactions['drug_name'].tolist() if i]
            # Look up names of interacting drugs from ChEMBL; add to previous list.
            interactions2 = self.chembl_interactions[self.chembl_interactions['hgnc_symbol'] == target['HGNC Name']]
            interactions_names.extend([str(i).upper() for i in interactions2['MOLECULE_NAME'].tolist() if i])
            # Check all the known drugs with warnings against the interacting list.
            if list(set(drugs_with_warnings_names).intersection(interactions_names)):
                return True
        
        def mouse_phenotype_association_severity():
            def contain(s, subs):
                # Does string s contain any of substrings subs?)
                return any(i in s for i in subs)
            
            # Do any of the jaxlab phenotype descriptors for this target contain any of the substrings identified as 
            # pointing to severe phenoptyes we want to give warnings about?
            jaxlab = self.jackson_lab_data[self.jackson_lab_data['HGNC Name'] == target['HGNC Name']]
            if any(contain(i, self.mouse_severe_phenotypes_patterns) for i in jaxlab['Term'].tolist()):
                return True
        
        def cancer_driver():
            if target['Mutational cancer driver genes'] or target['COSMIC somatic mutations in cancer genes']:
                return True
        
        def essential_gene():
            if target['Core fitness gene'] or target['CRISPR-screened core fitness gene'] \
                    or target['OGEE human essential gene'] in ['Essential', 'Conditional']:
                return True
        
        def single_gene_disorder():
            # Query the rare disease associations data for disease associations from OrphaNet (which are single-gene 
            # rare diseases)
            df = self.monogenic_disease_association_data[
                self.monogenic_disease_association_data['HGNC Name'] == target['HGNC Name']]
            if len(df.index):
                return True
        
        def hpa_enrichment_in_high_risk_tissues():
            # If this target is tissue enriched, group enriched or tissue enhanced in any HPA tissue marked by the user 
            # as a high risk site for overexpression 
            if [i for i in hpa_all_tissues if
                target["HPA " + i] and i in list(self.disease_profile['hpa']['off_target_high_risk'])]:
                return True
        
        def hpa_enrichment_in_low_risk_tissues():
            # If this target is tissue enriched, group enriched or tissue enhanced in any HPA tissue marked by the user 
            # as a low risk site for overexpression 
            if [i for i in hpa_all_tissues if
                target["HPA " + i] and i in list(self.disease_profile['hpa']['off_target_low_risk'])]:
                return True
        
        def hpo_associations_in_high_risk_tissues():
            # If this target is a gene marked as being associated with abnormalities in any HPO tissue marked by the 
            # user as an unacceptable site for abnormalities, return a 1 point result here. 
            if [i for i in hpo_annotation_names if
                target[i] and i in list(self.disease_profile['hpo']['high_risk'])]:
                return True
        
        def has_phase4_drug():
            drugs = self.existing_drug_data[self.existing_drug_data['HGNC Name'] == target['HGNC Name']]
            if 'Phase IV' in drugs['Max clinical phase'].tolist():
                return True
        
        def has_phase3_drug():
            drugs = self.existing_drug_data[self.existing_drug_data['HGNC Name'] == target['HGNC Name']]
            if 'Phase III' in drugs['Max clinical phase'].tolist():
                return True
        
        def has_existing_drug():
            drugs = self.existing_drug_data[self.existing_drug_data['HGNC Name'] == target['HGNC Name']]
            if "No existing drugs found" not in drugs['Drug name']:
                return True
        
        if has_existing_drug():
            if not boxed_warnings():
                if has_phase4_drug():
                    return 1
                if has_phase3_drug():
                    return 2
        if has_insufficient_info():
            return 6
        
        points = 0
        if on_target_boxed_warnings():
            points += 1
        if mouse_phenotype_association_severity():
            points += 1
        if cancer_driver():
            points += 2
        if essential_gene():
            points += 1
        if single_gene_disorder():
            points += 1
        if hpa_enrichment_in_high_risk_tissues():
            points += 2
        if hpa_enrichment_in_low_risk_tissues():
            points += 1
        if hpo_associations_in_high_risk_tissues():
            points += 1
        
        if on_target_withdrawn_drug():
            return 5
        
        if points >= 2:
            return 4
        elif points > 0:
            return 3
        else:
            return 2
    
    def get_homologs(self, uniprot_id):
        # Technically, paralog is the more correct term here, but even that refers to a qualitative rather than the
        # quantitative relationship we're looking at here. 
        # Shall we look up genes in the same family? That would be a good start, but we can do better - and use a 
        # process that falls more under our control and understanding.
        # I've made a dataset of identity scores from BLASTp alignments of all genes in the human proteome against
        # each other. The identity scores and e-values are in self.identity_scores and self.evalues respectively.
        threshold = self.api_keys['seq_similarity_threshold']
        if uniprot_id:
            if uniprot_id in self.identity_scores:
                idents = self.identity_scores[uniprot_id]
                idents_over_threshold = [i for i in idents if idents[i] >= threshold]
                # Remember - this figure is not symmetrical. I'll need to get it the other way too. 
                # Simplest thing would be to select only genes where the reciprocating value is ALSO above the 
                # threshold. This is a tad conservative, but not excessively so - and it's straightforward to code.
                # Oh - also, remember that this is basically a sparse matrix represented as a dict of dicts, so 
                # I'll have to use get() first, otherwise I'll get exceptions right quick. 
                inverse_idents = [i for i in idents_over_threshold if self.identity_scores.get(i).get(uniprot_id)]
                homologs = [i for i in inverse_idents if self.identity_scores[i][uniprot_id] >= threshold]
                # That will include the target uniprot ID; remove it.
                while uniprot_id in homologs:
                    homologs.remove(uniprot_id)
                return homologs
            else:
                return []
        else:
            return []
    
    def druggability_bucket(self, target):
        # Score for small-molecule druggability.
        
        def is_protein():
            if target['RNA class']:
                if "protein_coding" in target['RNA class'].split(", "):
                    return True
        
        def is_tclin_or_tchem(uniprot_id):
            # According to the Pharos website, there are also cases where targets are annotated as being Tclin/Tchem
            # (and therefore high priority targets for drug development) by manual curation. Include those here.
            if uniprot_id in self.pharos_ligand_data_library:
                pharos_data = self.pharos_ligand_data_library[uniprot_id]
                if pharos_data['Pharos'] in ["Tclin", "Tchem"]:
                    return True
        
        def has_ligand(uniprot_id):
            # Does this target have a ligand of any kind?
            # We may ask this question of genes that are not in our supplied list of targets, though that should be the
            # first place we look. It's plausible that the same homolog may come up multiple times too
            # We now keep a dict full of the previous results of queries to Pharos, to get around/streamline in those
            # situations.
            if uniprot_id in self.pharos_ligand_data_library:
                if self.pharos_ligand_data_library.get(uniprot_id).get('ChEMBL ligand') or \
                        self.pharos_ligand_data_library.get(uniprot_id).get('ChEMBL drug') or \
                        self.pharos_ligand_data_library.get(uniprot_id).get('ChEMBL low-activity ligand'):
                    return True
            elif target['Ligand'] or target['Endogenous ligand']:
                return True
            else:
                pharos_data = self.download_pharos_drug_and_ligand_count(uniprot_id)
                if pharos_data:
                    if pharos_data.get('ChEMBL ligand') or pharos_data.get('ChEMBL low-activity ligand'):
                        return True
        
        def has_high_activity_ligand(uniprot_id):
            # Similar to has_ligand, but restricts to ligands passing activity thresholds only.
            if uniprot_id in self.pharos_ligand_data_library:
                pharos_data = self.pharos_ligand_data_library[uniprot_id]
            else:
                pharos_data = self.download_pharos_drug_and_ligand_count(uniprot_id)
            if pharos_data:
                if pharos_data.get('ChEMBL ligand'):
                    return True
        
        def has_druggable_protein_class(uniprot_id):
            # Similar to has_ligand, but looks at whether the target has a druggable class listed in Pharos.
            if uniprot_id in self.pharos_ligand_data_library:
                if self.pharos_ligand_data_library.get(uniprot_id).get('ChEMBL ligand'):
                    return True
            else:
                pharos_data = self.download_pharos_drug_and_ligand_count(uniprot_id)
                if pharos_data:
                    if pharos_data.get('Druggable class'):
                        return True
        
        def has_endogenous_ligand():
            if target['Endogenous ligand']:
                return True
        
        def has_small_molecule_ligand():
            # We just checked if a target has any of a wide variety of ligands, but we only want to give higher scores 
            # in cases where there are non-endogenous ligands. This is easy to check when I have the right data: 
            # simply check that we have more total ligands than exclusively endogenous ligands.
            if target['Ligand'] > target['Endogenous ligand']:
                return True
        
        def homolog_has_ligand(uniprot_id):
            for hlog in self.get_homologs(uniprot_id):
                if has_ligand(hlog):
                    return True
        
        def homolog_has_high_activity_ligand(uniprot_id):
            for hlog in self.get_homologs(uniprot_id):
                if has_high_activity_ligand(hlog):
                    return True
        
        def has_structure(uniprot_id):
            if self.cansar_data[self.cansar_data['uniprot accession'] == uniprot_id]['number of 3d structure'].any():
                return True
        
        def has_druggable_pocket(uniprot_id):
            # if self.targets[self.targets['HGNC Name'] == name]['number of 3d structure druggable'].values[0]:
            if uniprot_id:
                if self.cansar_data[self.cansar_data['uniprot accession'] == uniprot_id][
                    'number of 3d structure druggable'].any():
                    return True
                # Should that fail to find anything, we can also check EBI's DrugEBIlity. 
                # drugebility = self.query_drugebility(uniprot_id)
                drugebility = self.query_ebi_drugebility_data(uniprot_id)
                if len(drugebility):
                    # if drugebility['EBI Druggable'] >= 0.8:
                    if len(drugebility[drugebility['DRUGGABLE'] >= 0.8]):
                        return True
        
        def homolog_has_structure(uniprot_id):
            for hlog in self.get_homologs(uniprot_id):
                if has_structure(hlog):
                    return True
        
        def homolog_has_druggable_pocket(uniprot_id):
            for hlog in self.get_homologs(uniprot_id):
                if has_druggable_pocket(hlog):
                    return True
        
        def family_has_structure(gene_family_id):
            for f in self.split_gene_families(gene_family_id):
                if self.gene_family_data[f].get('Approved Symbol'):  # All data outside of my control is just terrible
                    family_genes = self.gene_family_data[f]['Approved Symbol'].tolist()
                    if self.cansar_data[self.cansar_data['gene name'].isin(family_genes)][
                        'number of 3d structure'].any():
                        return True
        
        def family_has_druggable_pocket(gene_family_id):
            for f in self.split_gene_families(gene_family_id):
                family_genes = self.gene_family_data[f]['Approved Symbol'].tolist()
                if self.cansar_data[self.cansar_data['gene name'].isin(family_genes)][
                    'number of 3d structure druggable'].any():
                    return True
        
        def family_has_ligand(gene_family_id):
            for f in self.split_gene_families(gene_family_id):
                if self.gene_family_data[f].get('Approved Symbol'):  # All data outside of my control is just terrible
                    family_genes = self.gene_family_data[f]['Approved Symbol'].tolist()
                    for genename in family_genes:
                        uid = self.get_uniprot_id_for_gene_name(genename)
                        if has_ligand(uid):
                            return True
        
        def family_has_high_activity_ligand(gene_family_id):
            for f in self.split_gene_families(gene_family_id):
                family_genes = self.gene_family_data[f]['Approved Symbol'].tolist()
                for genename in family_genes:
                    uid = self.get_uniprot_id_for_gene_name(genename)
                    if has_high_activity_ligand(uid):
                        return True
        
        # Now we can stack them all up into our decision tree for this bucket. 
        name = target['HGNC Name']
        uniprot_id = target['Uniprot ID']
        gene_family_id = target['Gene family ID']
        if name and uniprot_id:
            if is_tclin_or_tchem(uniprot_id):
                return 1
            if has_ligand(uniprot_id):
                if has_small_molecule_ligand():
                    if has_high_activity_ligand(uniprot_id):
                        return 1
                    elif has_endogenous_ligand():
                        return 8
                    else:
                        return 5
                else:
                    return 8
            elif homolog_has_ligand(uniprot_id):
                if homolog_has_high_activity_ligand(uniprot_id):
                    return 2
                else:
                    return 6
            if has_structure(uniprot_id):
                if has_druggable_pocket(uniprot_id):
                    return 3
            if homolog_has_structure(uniprot_id) and homolog_has_druggable_pocket(uniprot_id):
                return 4
            if gene_family_id:
                if family_has_structure(gene_family_id):
                    if family_has_druggable_pocket(gene_family_id):
                        return 7
                if family_has_ligand(gene_family_id):
                    if family_has_high_activity_ligand(gene_family_id):
                        return 9
                    else:
                        return 10
            if has_druggable_protein_class(uniprot_id):
                return 11
            if has_structure(uniprot_id):
                return 12
            if is_protein():
                return 13
        # If nothing else has come through up to this point... 
        if is_protein():
            return 13
        else:
            return 14
    
    def feasibility_bucket(self, target):
        # Similar to other bucketing functions, I'll make a set of sub-functions that make individual component calls
        # and figure out a final output from there.
        def has_assays():
            if target['Assays']:
                return True
        
        def has_antibodies():
            # Pharos help indicates >=50 is a reasonable cutoff here.
            # We don't get an actual number here though, it seems, so this should be presence/absence too.
            if target['Antibodies']:
                return True
        
        def has_protein_structure_model():
            if target['Uniprot ID']:
                if len(self.pdb_data[self.pdb_data['Uniprot ID'] == target['Uniprot ID']]):
                    return True
        
        def has_protein_protein_interaction_data():
            if target['HGNC Name']:
                if len(self.ppi_data[self.ppi_data['HGNC Name'] == target['HGNC Name']]):
                    return True
        
        def has_protein_drug_interaction_data():
            if target['HGNC Name']:
                dgidb = self.dgidb_interactions[self.dgidb_interactions['gene_name'] == target['HGNC Name']]
                chembl = self.chembl_interactions[self.chembl_interactions['hgnc_symbol'] == target['HGNC Name']]
                if len(dgidb) > 0 or len(chembl) > 0:
                    return True
        
        def has_literature():
            if target.get("Literature"):
                if target['Literature'] > 0:
                    return True
            if target.get("Pharos literature"):
                if target['Pharos literature'] > 0:
                    return True
        
        def has_data_of_type(datatype, threshold=0):
            if target.get(datatype):
                if target[datatype] > threshold:
                    return True
        
        score = 0
        if has_assays():
            score += 1
        if has_antibodies():
            score += 1
        if has_protein_structure_model():
            score += 1
        if has_protein_protein_interaction_data():
            score += 1
        if has_literature():
            score += 1
        # if has_protein_drug_interaction_data:
        #     score += 1
        for t in self.feasibility_decision_data_types:
            if has_data_of_type(t):
                score += 1
        # We want to make it so that no information = 0, but all information = 1.  
        if score:
            return 8 - score
        else:
            return 0
    
    def antibodyability_bucket(self, target):
        # This is a little more straightforward than other bucketing functions, since everything we're deciding on is
        # already present (able to be gathered into a single sheet, no less) and requires very little further number
        # crunching.
        
        def secreted_strong_evidence():
            # This appears to be based on several different means of identifying tags for secreted proteins from the
            # HPA, so we can consider it a reasonably strong indicator. See
            # https://www.proteinatlas.org/humanproteome/secretome 
            if target['VerSeDa secretome membership']:
                return True
            # Also check if the HPA contains predictions of secreted proteins - if there is more than one method 
            # predicting secretion, this can be considered stronger evidence.
            if sum([target[i] for i in self.hpa_predicted_secreted_protein_subclasses if target[i]]) >= 2:
                return True
            # Also check if it's on record as localising to the Golgi apparatus with sufficient evidence.
            if search("Golgi apparatus", str(target['Subcellular location summary'])) \
                    and target['Subcellular location verification'] in ['approved', 'enhanced']:
                return True
        
        def ecm_strong_evidence():
            if target['Extracellular matrix component']:
                return True
        
        def membrane_strong_evidence():
            if target['Surfaceome membership (human)'] and target[
                'Surfaceome membership confidence (human)'] == "1 - high confidence":
                return True
            # Also check if the HPA contains predictions of membrane proteins - if there is more than one method 
            # predicting membrane location, this can be considered stronger evidence.
            if sum([target[i] for i in self.hpa_predicted_membrane_protein_subclasses if target[i]]) >= 2:
                return True
            # Also check if it's on record as localising to the Golgi apparatus with sufficient evidence.
            if search("plasma membrane", str(target['Subcellular location summary'])) \
                    and target['Subcellular location verification'] in ['approved', 'enhanced']:
                return True
        
        def secreted_ecm_membrane_weak_evidence():
            # Surfaceome membership (human): already checked for stronger evidence - this leaves only weaker evidence.
            # Surfaceome membership (mouse): this as a whole can be taken as weaker evidence.
            # Subcellular location summary: already checked for stronger evidence
            if target['Surfaceome membership (human)'] or \
                    target['Surfaceome membership (mouse)'] or \
                    target['Predicted membrane protein'] or \
                    target['Predicted secreted protein'] or \
                    [target[i] for i in self.hpa_predicted_membrane_protein_subclasses if target[i]] or \
                    [target[i] for i in self.hpa_predicted_secreted_protein_subclasses if target[i]] or \
                    search("plasma membrane|Golgi apparatus", str(target['Subcellular location summary'])):
                # Add further conditions as data become available...
                return True
        
        def cytoplasm():
            if target['Main subcellular locations']:
                if search("cytosol", str(target['Main subcellular locations'])):
                    return True
        
        def intracellular_compartment():
            # Basically any value in this field is now acceptable here
            if target['Main subcellular locations']:
                return True
        
        if secreted_strong_evidence():
            return 1
        elif ecm_strong_evidence():
            return 2
        elif membrane_strong_evidence():
            return 3
        elif secreted_ecm_membrane_weak_evidence():
            return 4
        elif cytoplasm():
            return 5
        elif intracellular_compartment():
            return 6
        else:
            return 7
    
    def new_modality_bucket(self, target, druggability, antibodyability):
        # Is this target a non-coding RNA? If so, it's the top priority for these modalities.
        # For a full list of RNA classes found in this field, see https://ensembl.org/Help/Faq?id=468
        if target['RNA class']:
            if "protein_coding" not in target['RNA class'].split(", "):
                return 1
        if 'New modality check' in target:
            if self.target['New modality check']:
                # If the user has supplied a value here, take it as confirmation that we should do this check.
                if druggability > 8 and antibodyability > 4:
                    return 1
                elif druggability > 8 or antibodyability > 4:
                    return 2
                else:
                    return 3
        else:
            return 4
    
    def tissue_engagement_bucket(self, target):
        def is_housekeeping_gene():
            if target['Housekeeping gene']:
                return True
        
        def enriched_in_high_risk_tissue():
            if [i for i in hpa_all_tissues if
                target["HPA " + i] and i in list(self.disease_profile['hpa']['off_target_high_risk'])]:
                return True
        
        def enriched_in_low_risk_tissue():
            if [i for i in hpa_all_tissues if
                target["HPA " + i] and i in list(self.disease_profile['hpa']['off_target_low_risk'])]:
                return True
        
        def not_enriched_in_any_tissue():
            if not [i for i in hpa_all_tissues if target["HPA " + i]]:
                return True
        
        def enriched_in_target_tissue():
            if [i for i in hpa_all_tissues if
                target["HPA " + i] and i in list(self.disease_profile['hpa']['on_target'])]:
                return True
        
        if enriched_in_high_risk_tissue():
            return 5
        if enriched_in_low_risk_tissue():
            return 4
        if is_housekeeping_gene():
            return 3
        if not_enriched_in_any_tissue():
            return 2
        if enriched_in_target_tissue():
            return 1


class PandasModel(QtCore.QAbstractTableModel):
    """
    Class to populate a table view with a pandas dataframe
    Modified from: http://stackoverflow.com/a/31557937
    """
    
    def __init__(self, data, parent=None, checkbox_col=False):
        QtCore.QAbstractTableModel.__init__(self, parent)
        data = self.data_type_cleanup(data)
        self._data = data.copy()
        if checkbox_col:
            self._data['Checked'] = False
            # data["Checked"] = QVariant(Qt.Unchecked)
            # Make that the first column
            self._data = self._data[['Checked'] + [i for i in data.columns.tolist() if i != "Checked"]]
            # print(data.columns.tolist())
        # Because we have copied the data, we can make a few modifications to these data frames in order to, e.g., 
        # make stuff look better, do sorts...
        # Replace None values with an empty string
        # self._data = self.replace_none_values(self._data)
        self._data = self._data.fillna('')
        self.orig_data = self._data.copy()
        self.headers = data.columns.values
        # self.print_dtypes()
    
    def data_type_cleanup(self, data):
        # Try to set some data types based on column values
        data = data.copy()
        try:
            for col in data.columns.tolist():
                if col in bool_cols:
                    # data[col] = data[col].astype(bool) 
                    # data[col] = data[col].map({1: "True", 0: None})
                    data[col] = ["True" if i else None for i in data[col].tolist()]
                else:
                    # Does this not match non-numeric characters? If so, 0 characters may not be recognised by the sort.
                    # data[col] = data[col].map({'0': 0})  # Removes EVERYTHING
                    data[col] = to_numeric(data[col], errors='ignore')
        except Exception as e:
            print('Caught exception in worker thread at model cleanup:')
            print_exc()
        return data
        
        # def flags(self, index):
    
    #     if 'Checked' in self._data.columns.tolist() and index.column() == 0:
    #         return Qt.ItemIsUserCheckable
    #     else:
    #         return Qt.ItemIsEnabled
    
    def replace_none_values(self, df):
        mask = df.applymap(lambda x: x is None or nan or 'nan')
        cols = df.columns[mask.any()]
        for col in df[cols]:
            df.loc[mask[col], col] = ''
        return df
    
    def rowCount(self, parent=None):
        return len(self._data.values)
    
    def columnCount(self, parent=None):
        return self._data.columns.size
    
    def data(self, index, role=QtCore.Qt.DisplayRole):
        try:
            if index.isValid():
                # if self._data.columns.tolist()[index.column()] == "Checked":
                #     if role == Qt.CheckStateRole:
                #         if self.dataObj[index.row()] == 0:
                #             return QVariant(Qt.Unchecked)
                #     else:
                #         return QVariant(Qt.Checked)
                if role == Qt.CheckStateRole and self._data.columns.tolist()[index.column()] == "Checked":
                    # return QVariant(Qt.Checked)
                    # return self._data.ix[index.row(), 'Checked']
                    if self._data.ix[index.row(), 'Checked']:
                        return QVariant(Qt.Checked)
                    else:
                        return QVariant(Qt.Unchecked)
                elif role == QtCore.Qt.DisplayRole:
                    return str(self._data.values[index.row()][index.column()])
            return None
        except Exception as e:
            print('Caught exception in worker thread:')
            print(e)
            print("Here's the data. That might be helpful.")
            print(self._data)
            print_exc()
    
    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self._data.columns[col]
        return None
    
    def sort(self, col, order):
        dtypes = {}
        
        print(self._data[self.headers[col - 1]].head())
        print(self._data[self.headers[col - 1]].dtypes)
        
        self.layoutAboutToBeChanged.emit()
        # Yes, this is technically the wrong way round. We want descending order on first click. 
        # Try converting to numeric - if there is any issue, just sort as strings.
        try:
            sort_by = self.headers[col - 1]
            print("Sorting by " + sort_by)
            self._data = self._data.sort_values(by=sort_by, ascending=False)
            if order == 1:
                self._data = self._data.sort_values(by=sort_by, ascending=True)
            self.layoutChanged.emit()
        except Exception as e:
            print('Caught exception in worker thread at model sort:')
            print_exc()
    
    def restrict(self, seriesnums):
        self.layoutAboutToBeChanged.emit()
        self._data = self.orig_data.copy()  # We do this to prevent multiple searches producing unexpected results
        self._data = self._data.loc[self._data['series'].isin(seriesnums)].reset_index(drop=True)
        self.layoutChanged.emit()
    
    def clear(self):
        self.layoutAboutToBeChanged.emit()
        self._data = self.orig_data.copy()
        self.layoutChanged.emit()
    
    def delete_everything(self):
        self.layoutAboutToBeChanged.emit()
        self._data.drop(self._data.index, inplace=True)
        self.orig_data.drop(self.orig_data.index, inplace=True)
        # When you clear everything out, you're still left with column headers. No idea how to remove them, and next to
        # no useful documentation is readily available. 
        # for i in range(len(self.headers)):
        #     self.removeColumn(0)
        self.layoutChanged.emit()
    
    def get_matching_indices(self, values, col="series"):
        # Get row indexes of data matching some set indices
        icol = self.get_column_number_for_column(col)
        indexes = [self.index(r, icol) for r in self._data[self._data[col].isin(values)].index.tolist()]
        return indexes
        # return self._data[self._data[col].isin(values)].index
    
    def get_matching_indices_in_original_data(self, values, col="series"):
        # Get row indexes of data matching some set indices
        icol = self.get_column_number_for_column(col)
        indexes = [self.index(r, icol) for r in self.orig_data[self.orig_data[col].isin(values)].index.tolist()]
        return indexes
    
    def get_column_number_for_column(self, col="series"):
        # Sometimes we want the numeric index of a column
        return self._data.columns.tolist().index(col)
    
    def check_rows(self, series_nos, check=True):
        # Also, uncheck rows.
        try:
            for ind in self.get_matching_indices(series_nos):
                self._data.at[ind.row(), 'Checked'] = check
                self.dataChanged.emit(ind, ind)
            # Should also set the original data, in case some restricting were to happen
            for ind in self.get_matching_indices_in_original_data(series_nos):
                self.orig_data.at[ind.row(), 'Checked'] = check
        except Exception as e:
            print('Caught exception in worker thread:')
            print_exc()
    
    def num_currently_checked(self):
        if "Checked" in self._data.columns.tolist():
            return len(self._data[self._data['Checked'] == True].index)
        else:
            return 0
    
    def checked_series_nums(self):
        if "Checked" in self._data.columns.tolist():
            return self._data[self._data['Checked'] == True]['series'].tolist()
        else:
            return 0
    
    def all_series_nums(self):
        return self._data['series'].tolist()


class Ui_ViewDataDialog(object):
    def setupUi(self, ViewDataDialog):
        try:
            ViewDataDialog.setObjectName("ViewDataDialog")
            ViewDataDialog.resize(1200, 800)
            # ViewDataDialog is a reference to the parent class, which lets us get at the excel sheet I loaded in there
            # self.sheetnames = list(ViewDataDialog.data.keys())
            self.sheetnames = [i for i in sheet_order if i not in unacceptably_large_sheets and i != "Bucket guide"]
            
            self.verticalLayout_4 = QtWidgets.QVBoxLayout(ViewDataDialog)
            self.verticalLayout_4.setObjectName("verticalLayout_4")
            self.label = QtWidgets.QLabel(ViewDataDialog)
            self.label.setObjectName("label")
            self.verticalLayout_4.addWidget(self.label)
            self.horizontalLayout = QtWidgets.QHBoxLayout()
            self.horizontalLayout.setObjectName("horizontalLayout")
            self.lineEdit = QtWidgets.QLineEdit(ViewDataDialog)
            self.lineEdit.setObjectName("lineEdit")
            self.horizontalLayout.addWidget(self.lineEdit)
            self.comboBox_searchType = QtWidgets.QComboBox(ViewDataDialog)
            self.comboBox_searchType.setMinimumSize(QtCore.QSize(200, 0))
            self.comboBox_searchType.setObjectName("comboBox_searchType")
            self.horizontalLayout.addWidget(self.comboBox_searchType)
            self.pushButton_search = QtWidgets.QPushButton(ViewDataDialog)
            self.pushButton_search.setMaximumSize(QtCore.QSize(51, 16777215))
            self.pushButton_search.setObjectName("pushButton_search")
            self.horizontalLayout.addWidget(self.pushButton_search)
            self.pushButton_clear = QtWidgets.QPushButton(ViewDataDialog)
            self.pushButton_clear.setMaximumSize(QtCore.QSize(41, 16777215))
            self.pushButton_clear.setObjectName("pushButton_clear")
            self.horizontalLayout.addWidget(self.pushButton_clear)
            self.verticalLayout_4.addLayout(self.horizontalLayout)
            self.tabWidget = QtWidgets.QTabWidget(ViewDataDialog)
            self.verticalLayout_4.addWidget(self.tabWidget)
            
            # This should give us a table view for every sheet in the excel file - with a few notable exceptions.
            self.tabs = dict()
            self.tableViews = dict()
            self.verticalLayouts = dict()
            # for sheet in [i for i in sheet_order if i not in unacceptably_large_sheets and i != "Bucket guide"]:
            #     if sheet in self.sheetnames:
            for sheet in self.sheetnames:
                print("Setting up tab for " + sheet)
                self.tabs[sheet] = QtWidgets.QWidget()
                self.tabs[sheet].setObjectName("tab")
                self.verticalLayouts[sheet] = QtWidgets.QVBoxLayout(self.tabs[sheet])
                self.verticalLayouts[sheet].setObjectName("verticalLayout")
                self.tableViews[sheet] = QTableView(self.tabs[sheet])
                self.tableViews[sheet].setObjectName("tableView")
                self.verticalLayouts[sheet].addWidget(self.tableViews[sheet])
                self.tableViews[sheet].raise_()
                self.tableViews[sheet].raise_()
                self.tabWidget.addTab(self.tabs[sheet], "")
                self.tableViews[sheet].setSortingEnabled(True)
                self.tableViews[sheet].setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
                self.tableViews[sheet].setModel(ViewDataDialog.models[sheet])
                # If a checkbox column is required for this tab, set it up here 
                # if 'Checked' in ViewDataDialog.models[sheet].orig_data.columns.tolist():
                # print("Trying to set delegate for " + sheet)
                # self.tableViews[sheet].openPersistentEditor()
                # self.tableViews[sheet].setItemDelegateForColumn(0, CheckBoxDelegate(self.tableViews[sheet]))
                # self.tableViews[sheet].setItemDelegateForColumn(0, OtherCheckboxDelegate(self.tableViews[sheet]))
                self.tableViews[sheet].resizeColumnsToContents()
            # This bit also wants to be kept separate, for ease of copy-pasting.
            search_types = ['HGNC Name', 'GeneID', 'Uniprot ID']
            self.comboBox_searchType.addItems(search_types)
            
            self.gridLayout = QtWidgets.QGridLayout()
            self.gridLayout.setObjectName("gridLayout")
            self.pushButton_checkAll = QtWidgets.QPushButton(ViewDataDialog)
            self.pushButton_checkAll.setMinimumSize(QtCore.QSize(61, 0))
            self.pushButton_checkAll.setMaximumSize(QtCore.QSize(61, 51))
            self.pushButton_checkAll.setObjectName("pushButton_checkAll")
            self.gridLayout.addWidget(self.pushButton_checkAll, 0, 0, 1, 1)
            self.verticalLayout = QtWidgets.QVBoxLayout()
            self.verticalLayout.setObjectName("verticalLayout")
            self.pushButton_checkSelected = QtWidgets.QPushButton(ViewDataDialog)
            self.pushButton_checkSelected.setEnabled(False)
            self.pushButton_checkSelected.setMinimumSize(QtCore.QSize(91, 0))
            self.pushButton_checkSelected.setMaximumSize(QtCore.QSize(91, 16777215))
            self.pushButton_checkSelected.setObjectName("pushButton_checkSelected")
            self.verticalLayout.addWidget(self.pushButton_checkSelected)
            self.pushButton_clearSelection = QtWidgets.QPushButton(ViewDataDialog)
            self.pushButton_clearSelection.setEnabled(False)
            self.pushButton_clearSelection.setMinimumSize(QtCore.QSize(91, 0))
            self.pushButton_clearSelection.setMaximumSize(QtCore.QSize(91, 16777215))
            self.pushButton_clearSelection.setObjectName("pushButton_clearSelection")
            self.verticalLayout.addWidget(self.pushButton_clearSelection)
            self.gridLayout.addLayout(self.verticalLayout, 0, 1, 1, 1)
            self.verticalLayout_2 = QtWidgets.QVBoxLayout()
            self.verticalLayout_2.setObjectName("verticalLayout_2")
            self.pushButton_unCheckSelected = QtWidgets.QPushButton(ViewDataDialog)
            self.pushButton_unCheckSelected.setEnabled(False)
            self.pushButton_unCheckSelected.setMinimumSize(QtCore.QSize(105, 0))
            self.pushButton_unCheckSelected.setMaximumSize(QtCore.QSize(105, 16777215))
            self.pushButton_unCheckSelected.setObjectName("pushButton_unCheckSelected")
            self.verticalLayout_2.addWidget(self.pushButton_unCheckSelected)
            self.pushButton_checkedOnly = QtWidgets.QPushButton(ViewDataDialog)
            self.pushButton_checkedOnly.setEnabled(False)
            self.pushButton_checkedOnly.setMinimumSize(QtCore.QSize(105, 0))
            self.pushButton_checkedOnly.setMaximumSize(QtCore.QSize(105, 16777215))
            self.pushButton_checkedOnly.setObjectName("pushButton_checkedOnly")
            self.verticalLayout_2.addWidget(self.pushButton_checkedOnly)
            self.gridLayout.addLayout(self.verticalLayout_2, 0, 2, 1, 1)
            spacerItem = QtWidgets.QSpacerItem(218, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
            self.gridLayout.addItem(spacerItem, 0, 3, 1, 1)
            self.pushButton_createReports = QtWidgets.QPushButton(ViewDataDialog)
            self.pushButton_createReports.setEnabled(False)
            self.pushButton_createReports.setMinimumSize(QtCore.QSize(91, 51))
            self.pushButton_createReports.setMaximumSize(QtCore.QSize(91, 51))
            self.pushButton_createReports.setObjectName("pushButton_createReports")
            self.gridLayout.addWidget(self.pushButton_createReports, 0, 4, 1, 1)
            self.verticalLayout_4.addLayout(self.gridLayout)
            
            self.retranslateUi(ViewDataDialog)
            self.tabWidget.setCurrentIndex(0)
            QtCore.QMetaObject.connectSlotsByName(ViewDataDialog)
        except Exception as e:
            print('Caught exception in worker thread:')
            print_exc()
    
    def retranslateUi(self, ViewDataDialog):
        _translate = QtCore.QCoreApplication.translate
        ViewDataDialog.setWindowTitle(_translate("ViewDataDialog", "Data viewer"))
        self.label.setText(_translate("ViewDataDialog", "Enter a list of genes to restrict data"))
        self.pushButton_search.setText(_translate("ViewDataDialog", "Search"))
        self.pushButton_clear.setText(_translate("ViewDataDialog", "Clear"))
        self.pushButton_checkAll.setText(_translate("ViewDataDialog", "Check\nall"))
        self.pushButton_checkSelected.setText(_translate("ViewDataDialog", "Check selected"))
        self.pushButton_clearSelection.setText(_translate("ViewDataDialog", "Clear selection"))
        self.pushButton_createReports.setText(_translate("ViewDataDialog", "Create reports"))
        self.pushButton_unCheckSelected.setText(_translate("ViewDataDialog", "Uncheck selected"))
        self.pushButton_checkedOnly.setText(_translate("ViewDataDialog", "Show checked only"))
        for sheet in self.sheetnames:
            self.tabWidget.setTabText(self.tabWidget.indexOf(self.tabs[sheet]), _translate("ViewDataDialog", sheet))


class DataWindow(QDialog):
    make_reports = pyqtSignal(list)
    
    def __init__(self, df=None, excel=None):
        super(DataWindow, self).__init__()
        # If given a list of data frames, feed that directly into the window viewer. 
        # If given a path to an excel file, read that into a list of data frames and feed that directly into the window
        # viewer.
        self.data = None
        self.models = dict()
        self.selected_series_nos = []
        if not df and not excel:
            w = QWidget()
            QMessageBox.critical(w, "Error", "Expected either list of dataframes or path to excel file!")
        elif df:
            self.data = df
            self.make_models()
        elif excel:
            self.read_excel_file(excel)
    
    def read_excel_file(self, excel):
        print("Reading excel file")
        self.data = read_excel(excel,
                               sheet_name=None,  # sheet_name=None gets us all sheets, in a dict
                               index_col=None,
                               converters={'gene_symbol': str,
                                           'chromosome_name': str})
        self.make_models()
    
    def make_models(self):
        for sheet in [i for i in self.data if i not in unacceptably_large_sheets and i != "Bucket guide"]:
            try:
                checkbox_col = False
                ls = self.data[sheet]['series'].tolist()
                if len(uniq(ls)) == len(ls):
                    checkbox_col = True
                self.models[sheet] = PandasModel(self.data[sheet], checkbox_col=checkbox_col)
            except Exception as e:
                print('Caught exception in worker thread:')
                print_exc()
        self.set_up_ui()
    
    def set_up_ui(self):
        self.ui = Ui_ViewDataDialog()
        self.ui.setupUi(self)
        self.setWindowIcon(QIcon('logo2.png'))
        self.ui.pushButton_search.clicked.connect(self.search_data)
        self.ui.pushButton_clear.clicked.connect(self.clear_search)
        self.ui.pushButton_checkAll.clicked.connect(self.check_all)
        self.ui.pushButton_clearSelection.clicked.connect(self.clear_selection)
        self.ui.pushButton_checkSelected.clicked.connect(self.check_selection)
        self.ui.pushButton_unCheckSelected.clicked.connect(self.uncheck_selection)
        # self.ui.pushButton_checkedOnly.clicked.connect(self.filter_to_checked_only)
        self.change_checkedOnly_button(True)
        self.ui.pushButton_createReports.clicked.connect(self.create_reports)
        # Set up 'clicked' (or whatever from the selection model) for the tableviews - selecting a row on one should 
        # select matching values across all.
        # try:
        #     for sheetname in self.ui.tableViews:
        #         tview = self.ui.tableViews[sheetname]
        #         tview.selectionModel().selectionChanged.connect(self.apply_selection_to_all_tables)
        #         # If that doesn't work, try this:
        #         # selectionModel = tview.selectionModel()
        #         # selectionModel.selectionChanged.connect(self.apply_selection_to_all_tables)
        # except Exception as e:
        #     print('Caught exception in worker thread:')
        #     print_exc()
        # Change of plan: 
        # When selection is changed, check if the 'clear selection' button should be activated or deactivated.
        try:
            self.ui.tabWidget.currentChanged.connect(self.activate_deactivate_clearSelection_button)
            for sheetname in self.ui.tableViews:
                tview = self.ui.tableViews[sheetname]
                tview.selectionModel().selectionChanged.connect(self.activate_deactivate_clearSelection_button)
                # If that doesn't work, try this:
                # selectionModel = tview.selectionModel()
                # selectionModel.selectionChanged.connect(self.apply_selection_to_all_tables)
        except Exception as e:
            print('Caught exception in worker thread:')
            print_exc()
    
    def activate_deactivate_clearSelection_button(self):
        # We're only doing this for the current tab - others are out of sight and out of mind. 
        # The 'check selected' and 'uncheck selected' buttons are along for the same ride, so do those at the same time.
        current_tab = self.ui.tabWidget.tabText(self.ui.tabWidget.currentIndex())
        tview = self.ui.tableViews[current_tab]
        inds = tview.selectedIndexes()
        num_checked = self.models[current_tab].num_currently_checked()
        if inds:
            self.ui.pushButton_clearSelection.setEnabled(True)
            self.ui.pushButton_checkSelected.setEnabled(True)
        else:
            self.ui.pushButton_clearSelection.setEnabled(False)
            self.ui.pushButton_checkSelected.setEnabled(False)
        if num_checked > 0:
            self.ui.pushButton_checkedOnly.setEnabled(True)
            self.ui.pushButton_createReports.setEnabled(True)
        else:
            self.ui.pushButton_checkedOnly.setEnabled(False)
            self.ui.pushButton_createReports.setEnabled(False)
        if inds or num_checked > 0:
            self.ui.pushButton_unCheckSelected.setEnabled(True)
        else:
            self.ui.pushButton_unCheckSelected.setEnabled(False)
    
    def search_data(self):
        # Emit a signal containing the parsed list; the table models can pick it up and act appropriately. 
        names = split("[\s\n:;,]+", str(self.ui.lineEdit.text()))
        if len(names) >= 1:
            col = str(self.ui.comboBox_searchType.currentText())
            seriesnums = self.data['Basic information'][self.data['Basic information'][col].isin(names)]['series']
            for sheet in self.data:
                self.models[sheet].restrict(seriesnums)
    
    def clear_search(self):
        self.ui.lineEdit.setText("")
        for sheet in self.data:
            self.models[sheet].clear()
        # Got to reset the checkedOnly button when this is clicked, 'cause it invokes this.
        self.change_checkedOnly_button(True)
    
    def apply_selection_to_all_tables(self, selected, deselected):
        try:
            # items = selected.selectedItems()
            # print([str(i) for i in items])
            current_tab = self.ui.tabWidget.tabText(self.ui.tabWidget.currentIndex())
            # print("I think a selection was changed on tab " + str(current_tab))
            # pprint(selected)
            # print("Columns:")
            # print([i.column() for i in selected.indexes()])
            # print("All values:")
            # print([i.data() for i in selected.indexes()])
            # print("Filtered series values:")
            # print([i.data() for i in selected.indexes() if i.column() == 0])
            
            # Get the data model for the current tab - we can use it to unambiguously look up the index of the column 
            # with the name 'series'
            series_index = self.models[current_tab].get_column_number_for_column()
            series_nos_selected = [int(i.data()) for i in selected.indexes() if i.column() == series_index]
            series_nos_deselected = [int(i.data()) for i in deselected.indexes() if i.column() == series_index]
            self.selected_series_nos = [int(i) for i in
                                        list(set(self.selected_series_nos).union(set(series_nos_selected)))]
            self.selected_series_nos = sorted(
                [int(i) for i in list(set(self.selected_series_nos) - set(series_nos_deselected))])
            # print("Better filtered series values:")
            # print(series_nos_selected)
            # print(series_nos_deselected)
            # Now loop through all other tables and select (and deselect) rows where the series value matches the ones
            # we just pulled out.
            for sheetname in [i for i in list(self.data.keys()) if i is not current_tab]:
                df = self.data[sheetname]
                # selected_indexes = [self.models[sheetname].index(r, 0) for r in 
                #                     df[df['series'].isin(series_nos_selected)].index.tolist()]
                # self.selected_series_nos]
                selected_indexes = self.models[sheetname].get_matching_indices(series_nos_selected)
                # print("Need to select these indices in sheet " + sheetname + ":")
                # print(selected_indexes)
                # self.ui.tableViews[sheetname].selectionModel().select(selected_indexes)
                mode = QItemSelectionModel.Select | QItemSelectionModel.Rows
                for ind in selected_indexes:
                    self.ui.tableViews[sheetname].selectionModel().select(ind, mode)
                    # self.ui.tableViews[sheetname].selectionModel().setSelected(i, mode)
                # deselected_indexes df[self.data[sheetname]['series'] in series_nos_deselected].index.tolist()
                # self.ui.tableViews[sheetname].selectionModel().Deselect(deselected_indexes)
                mode = QItemSelectionModel.Deselect | QItemSelectionModel.Rows
                # deselected_indexes = [self.models[sheetname].index(r, 0) for r in
                #                       df[df['series'].isin(series_nos_deselected)].index.tolist()]
                deselected_indexes = self.models[sheetname].get_matching_indices(series_nos_deselected)
                for ind in deselected_indexes:
                    self.ui.tableViews[sheetname].selectionModel().select(ind, mode)
        except Exception as e:
            print('Caught exception in worker thread:')
            print_exc()
    
    def check_all(self):
        # Check everything that's currently on the screen
        sheet = 'Basic information'
        selected_series = self.models[sheet].all_series_nums()
        print("CHECK ALL")
        print(selected_series)
        for sheetname in list(self.data.keys()):
            self.models[sheetname].check_rows(selected_series, check=True)
            self.ui.tableViews[sheetname].selectionModel().clear()
    
    def clear_selection(self):
        print("CLEAR SELECTION")
        # for sheetname in [i for i in list(self.data.keys())]:
        #     self.ui.tableViews[sheetname].selectionModel().clear()
        # Actually, apply it only to one tab at a time. 
        current_tab = self.ui.tabWidget.tabText(self.ui.tabWidget.currentIndex())
        self.ui.tableViews[current_tab].selectionModel().clear()
        self.activate_deactivate_clearSelection_button()
    
    def check_selection(self):
        current_tab = self.ui.tabWidget.tabText(self.ui.tabWidget.currentIndex())
        tview = self.ui.tableViews[current_tab]
        inds = tview.selectedIndexes()
        # Get the actual 'series' numbers for that. They need to be checked across all tables. There should be funcs 
        # built into the pandas model for this.
        series_index = self.models[current_tab].get_column_number_for_column()
        selected_series = [int(i.data()) for i in inds if i.column() == series_index]
        for sheetname in list(self.data.keys()):
            self.models[sheetname].check_rows(selected_series, check=True)
        self.ui.tableViews[current_tab].selectionModel().clear()
    
    def uncheck_selection(self):
        current_tab = self.ui.tabWidget.tabText(self.ui.tabWidget.currentIndex())
        tview = self.ui.tableViews[current_tab]
        inds = tview.selectedIndexes()
        # Get the actual 'series' numbers for that. They need to be checked across all tables. There should be funcs 
        # built into the pandas model for this.
        series_index = self.models[current_tab].get_column_number_for_column()
        selected_series = [int(i.data()) for i in inds if i.column() == series_index]
        for sheetname in list(self.data.keys()):
            self.models[sheetname].check_rows(selected_series, check=False)
        self.ui.tableViews[current_tab].selectionModel().clear()
    
    def filter_to_checked_only(self):
        sheet = 'Basic information'
        seriesnums = self.models[sheet].checked_series_nums()
        for sheet in self.data:
            self.models[sheet].restrict(seriesnums)
        self.change_checkedOnly_button(False)
    
    def change_checkedOnly_button(self, state):
        if state:
            self.ui.pushButton_checkedOnly.setText("Show checked only")
            self.ui.pushButton_checkedOnly.clicked.connect(self.filter_to_checked_only)
        else:
            self.ui.pushButton_checkedOnly.setText("Re-show all data")
            self.ui.pushButton_checkedOnly.clicked.connect(self.clear_search)
    
    def create_reports(self):
        # Get list of targets with checks
        sheet = 'Basic information'
        seriesnums = self.models[sheet].checked_series_nums()
        print("Selected these genes for reporting:")
        print(seriesnums)
        # Good! Now go back to the main thread and kick off the code that handles report generation. 
        self.make_reports.emit(seriesnums)


class TextInsertDialog(QDialog):
    def __init__(self, ):
        super(TextInsertDialog, self).__init__()
        self.ui = textinsert.Ui_Dialog()
        self.ui.setupUi(self)
        self.setWindowIcon(QIcon('logo2.png'))
        self.target_list = []
        self.ui.lineEdit.setText(str(len(self.target_list)) + " targets")
        self.ui.textEdit.textChanged.connect(self.parse_gene_list)
    
    def parse_gene_list(self):
        self.target_list = [i for i in split("[\s\n:;,]+", self.ui.textEdit.toPlainText()) if len(i) > 0]
        self.ui.lineEdit.setText(str(len(self.target_list)) + " targets")
    
    def get_values(self):
        return self.target_list, str(self.ui.buttonGroup.checkedButton().text())


class DiseaseProfileSelectorDialog(QDialog):
    def __init__(self, downloader, default_profile=None):
        # We need a reference to the downloader class so we can get access to the expression datasets to populate the 
        # fields in this dialog. 
        try:
            super(DiseaseProfileSelectorDialog, self).__init__()
            self.ui = disease_profile_dialog.Ui_Dialog()
            self.ui.setupUi(self)
            self.setWindowIcon(QIcon('logo2.png'))
            self.downloader = downloader
            
            self.ui.pushButton_offTarget_lowRisk_to_OnTarget_hpa.clicked.connect(
                lambda: self.move_between_lists(source_list=self.ui.listWidget_offTarget_lowRisk_hpa,
                                                destination_list=self.ui.listWidget_onTarget_hpa))
            self.ui.pushButton_onTarget_to_OffTarget_lowRisk_hpa.clicked.connect(
                lambda: self.move_between_lists(source_list=self.ui.listWidget_onTarget_hpa,
                                                destination_list=self.ui.listWidget_offTarget_lowRisk_hpa))
            
            self.ui.pushButton_offTarget_highRisk_to_offTarget_lowRisk_hpa.clicked.connect(
                lambda: self.move_between_lists(source_list=self.ui.listWidget_offTarget_highRisk_hpa,
                                                destination_list=self.ui.listWidget_offTarget_lowRisk_hpa))
            self.ui.pushButton_offTarget_lowRisk_to_offTarget_highRisk_hpa.clicked.connect(
                lambda: self.move_between_lists(source_list=self.ui.listWidget_offTarget_lowRisk_hpa,
                                                destination_list=self.ui.listWidget_offTarget_highRisk_hpa))
            
            self.ui.pushButton_highRisk_to_lowRisk_hpo.clicked.connect(
                lambda: self.move_between_lists(source_list=self.ui.listWidget_highRisk_hpo,
                                                destination_list=self.ui.listWidget_lowRisk_hpo))
            self.ui.pushButton_lowRisk_to_highRisk_hpo.clicked.connect(
                lambda: self.move_between_lists(source_list=self.ui.listWidget_lowRisk_hpo,
                                                destination_list=self.ui.listWidget_highRisk_hpo))
            self.populate_saved_profiles_list()
            self.ui.pushButton_reset_default.clicked.connect(self.set_default)
            # Now this stuff needs to be done in the right order, so that we open up on the previous setting.
            if default_profile:
                self.ui.comboBox_selectProfile.setCurrentIndex(
                    self.ui.comboBox_selectProfile.findText(default_profile['name']))
                self.existing_profile_selected(self.ui.comboBox_selectProfile.findText(default_profile['name']))
            else:
                self.set_default()
            self.ui.comboBox_selectProfile.currentIndexChanged.connect(self.existing_profile_selected)
            self.ui.lineEdit_newProfileName.textChanged.connect(self.activate_save_button)
            self.ui.pushButton_save.clicked.connect(self.save_profile)
        except Exception as e:
            print("Exception in disease profile thing:")
            print(e)
            print_exc()
    
    def set_default(self):
        # Put everything back in a starting position.
        self.ui.listWidget_onTarget_hpa.clear()
        self.ui.listWidget_offTarget_lowRisk_hpa.clear()
        self.ui.listWidget_offTarget_highRisk_hpa.clear()
        self.ui.listWidget_lowRisk_hpo.clear()
        self.ui.listWidget_highRisk_hpo.clear()
        for i in hpa_all_tissues:
            if i in ['Heart muscle', 'Liver', 'Kidney']:
                self.ui.listWidget_offTarget_highRisk_hpa.addItem(i)
            else:
                self.ui.listWidget_offTarget_lowRisk_hpa.addItem(i)
            self.ui.listWidget_offTarget_lowRisk_hpa.sortItems()
            self.ui.listWidget_offTarget_highRisk_hpa.sortItems()
        for i in hpo_annotation_names:
            if i in ['Cardiovascular', 'Hepatic', 'Renal']:
                self.ui.listWidget_highRisk_hpo.addItem(i)
            else:
                self.ui.listWidget_lowRisk_hpo.addItem(i)
            self.ui.listWidget_lowRisk_hpo.sortItems()
            self.ui.listWidget_highRisk_hpo.sortItems()
        self.ui.comboBox_selectProfile.setCurrentText("Default")
        self.ui.lineEdit_newProfileName.clear()
    
    def activate_save_button(self):
        # Is there anything written in the 'Save new profile' line edit? If so, activate the Save button. If not, 
        # deactivate it. This prevents users from clicking it in a circumstance that might cause trouble.
        try:
            if len(self.ui.lineEdit_newProfileName.text()) > 0:
                self.ui.pushButton_save.setEnabled(True)
            else:
                self.ui.pushButton_save.setEnabled(False)
        except Exception as e:
            print("Exception in disease profile thing:")
            print(e)
            print_exc()
    
    def populate_saved_profiles_list(self):
        try:
            profiles = [basename(i) for i in glob("disease_profiles/*.json")]
            if 'Default' not in profiles:
                profiles.append('Default')
            for i in profiles:
                profname = sub('.json', '', i)
                # profname = sub("_", " ", profname)
                self.ui.comboBox_selectProfile.addItem(profname)
        except Exception as e:
            print("Exception in disease profile thing:")
            print(e)
            print_exc()
    
    def make_some_json_output(self):
        try:
            return {
                "name": self.ui.lineEdit_newProfileName.text() or self.ui.comboBox_selectProfile.currentText(),
                "hpa": {
                    "on_target": self.list_items_in_listwidget(self.ui.listWidget_onTarget_hpa),
                    "off_target_low_risk": self.list_items_in_listwidget(self.ui.listWidget_offTarget_lowRisk_hpa),
                    "off_target_high_risk": self.list_items_in_listwidget(self.ui.listWidget_offTarget_highRisk_hpa)},
                "hpo": {
                    "low_risk": self.list_items_in_listwidget(self.ui.listWidget_lowRisk_hpo),
                    "high_risk": self.list_items_in_listwidget(self.ui.listWidget_highRisk_hpo)}
            }
        except Exception as e:
            print("Exception in disease profile thing:")
            print(e)
            print_exc()
    
    def save_profile(self):
        try:
            profile = self.make_some_json_output()
            with open("disease_profiles/" + profile['name'] + ".json", 'w') as outfile:
                dump(profile, outfile)
            # I should also update the combobox to have this name in it and selected!
            self.populate_saved_profiles_list()
            self.ui.lineEdit_newProfileName.clear()
            # self.ui.comboBox_selectProfile.setCurrentText(profile['name'])
            self.ui.comboBox_selectProfile.setCurrentIndex(self.ui.comboBox_selectProfile.findText(profile['name']))
        except Exception as e:
            print("Exception in disease profile thing:")
            print(e)
            print_exc()
    
    def existing_profile_selected(self, index):
        try:
            selected_profile = self.ui.comboBox_selectProfile.itemText(index)
            print("Selected " + selected_profile)
            if selected_profile == "Default":
                self.set_default()
            else:
                profdata = load(open("disease_profiles/" + selected_profile + ".json"))
                self.ui.listWidget_onTarget_hpa.clear()
                self.ui.listWidget_offTarget_lowRisk_hpa.clear()
                self.ui.listWidget_offTarget_highRisk_hpa.clear()
                self.ui.listWidget_lowRisk_hpo.clear()
                self.ui.listWidget_highRisk_hpo.clear()
                self.ui.listWidget_onTarget_hpa.addItems(profdata['hpa']['on_target'])
                self.ui.listWidget_offTarget_lowRisk_hpa.addItems(profdata['hpa']['off_target_low_risk'])
                self.ui.listWidget_offTarget_highRisk_hpa.addItems(profdata['hpa']['off_target_high_risk'])
                self.ui.listWidget_lowRisk_hpo.addItems(profdata['hpo']['low_risk'])
                self.ui.listWidget_highRisk_hpo.addItems(profdata['hpo']['high_risk'])
        except Exception as e:
            print("Exception in disease profile thing:")
            print(e)
            print_exc()
    
    def list_items_in_listwidget(self, listwidget):
        try:
            return [str(listwidget.item(i).text()) for i in range(listwidget.count())]
        except Exception as e:
            print("Exception in disease profile thing:")
            print(e)
            print_exc()
    
    def get_index_of_matching_value_in_listwidget(self, listwidget, value):
        try:
            c = 0
            for i in self.list_items_in_listwidget(listwidget):
                if i == value:
                    return c
                c += 1
            return None
        except Exception as e:
            print("Exception in disease profile thing:")
            print(e)
            print_exc()
    
    def moved_items(self, listwidget):
        try:
            move_these = [x.text() for x in listwidget.selectedItems()]
            leave_these = list(set(self.list_items_in_listwidget(listwidget)) - set(move_these))
            return move_these, leave_these
        except Exception as e:
            print("Exception in disease profile thing:")
            print(e)
            print_exc()
    
    def move_between_lists(self, source_list, destination_list):
        try:
            moved_items, remaining_items = self.moved_items(source_list)
            for i in moved_items:
                destination_list.addItem(i)
            source_list.clear()
            for i in remaining_items:
                source_list.addItem(i)
            source_list.sortItems()
            destination_list.sortItems()
        except Exception as e:
            print("Exception in disease profile thing:")
            print(e)
            print_exc()


class ReportGeneratorDialog(QDialog):
    # Two things here.
    # 1: Fire off the loop that goes and makes reports.
    # 2: Show progress in a dialog.
    
    # We WANT the main UI thread to be locked while doing this, so I won't be making this a separate thread. 
    #  However, I think I'll have to make the dialog have its own thread. 
    
    def __init__(self):
        try:
            # super(QThread, self).__init__()
            super(QDialog, self).__init__()
            self.ui = makereports_dialog.Ui_Dialog()
            self.ui.setupUi(self)
            self.setWindowIcon(QIcon('logo2.png'))
        except Exception as e:
            print("Exception in report gen function dialog thing:")
            print(e)
            print_exc()
    
    def update_label(self, gene_name):
        self.ui.label.setText("Generating report for target... " + gene_name)
    
    def update_progbar(self, pc):
        # Pass in a percentage figure
        pc = round(pc)  # This gives us an int, and also gets around any potential floating point-related strangeness. 
        self.ui.progressBar.setValue(pc)
        if pc == 100:
            self.ui.buttonBox.setEnabled(True)


class ReportGeneratorThreadThing(QThread):
    def __init__(self):
        super(QThread, self).__init__()
        self.dialog = ReportGeneratorDialog()
    
    def run(self):
        self.dialog.show()
        self.dialog.exec_()
    
    def update_label(self, gene_name):
        self.dialog.update_label(gene_name)
    
    def update_progbar(self, pc):
        self.dialog.update_progbar(pc)


class ReportGenerator:
    # gene = pyqtSignal(str)
    # progbar = pyqtSignal(int)
    
    def __init__(self, df, series, outfile):
        self.data = df
        self.series = series
        self.outfile = outfile
        # self.dialog_thread = ReportGeneratorDialog()
        # self.dialog_thread = ReportGeneratorThreadThing()
    
    def make_reports(self):
        try:
            # self.dialog_thread.start()
            # self.dialog_thread.show()
            i = 0
            for s in self.series:
                gene_name = \
                    self.data['Basic information'][self.data['Basic information']['series'] == s]['HGNC Name'].values[0]
                print(gene_name)
                # Write some temporary CSV files, invoke the report markdown thing, and then delete the temp files.
                # That's messy, but it's one of the quicker ways of doing this - especially when dealing with a lot of 
                # targets. 
                # I thought about doing some input (for the report thing) verification, but a better approach would be 
                # to make it robust against anything that's going to come out of here. Any issues should have been 
                # sorted out long ago by this point.
                temp_gtex_file = "temp/gtex.csv"
                temp_basic_info = "temp/basic.csv"
                temp_pharos = "temp/pharos.csv"
                temp_barres_mouse = "temp/barres_mouse.csv"
                temp_barres_human = "temp/barres_human.csv"
                temp_disease_assocs = "temp/disease_assocs.csv"
                temp_drugs = "temp/drugs.csv"
                temp_druggability = "temp/druggability.csv"
                temp_ab = "temp/ab.csv"
                temp_risk = "temp/risk.csv"
                temp_feas = "temp/feas.csv"
                
                self.data['GTEX'][self.data['GTEX']['series'] == s].to_csv(temp_gtex_file)
                self.data['Basic information'][self.data['Basic information']['series'] == s].to_csv(temp_basic_info)
                self.data['Pharos'][self.data['Pharos']['series'] == s].to_csv(temp_pharos)
                self.data['Barres mouse'][self.data['Barres mouse']['series'] == s].to_csv(temp_barres_mouse)
                self.data['Barres human'][self.data['Barres human']['series'] == s].to_csv(temp_barres_human)
                self.data['Common disease associations'][
                    self.data['Common disease associations']['series'] == s].to_csv(
                    temp_disease_assocs)
                self.data['Rare disease associations'][
                    self.data['Rare disease associations']['series'] == s].to_csv(
                    temp_disease_assocs)
                self.data['Existing drugs'][self.data['Existing drugs']['series'] == s].to_csv(temp_drugs)
                self.data['SM Druggability'][self.data['SM Druggability']['series'] == s].to_csv(temp_druggability)
                self.data['AB-ability'][self.data['AB-ability']['series'] == s].to_csv(temp_ab)
                self.data['Risk factors'][self.data['Risk factors']['series'] == s].to_csv(temp_risk)
                
                report_html_file = gene_name + "_report.html"
                weave("summary_report.md", doctype='md2html', informat='markdown', output=report_html_file)
                
                for f in [temp_gtex_file, temp_basic_info, temp_pharos, temp_barres_mouse, temp_barres_human,
                          temp_disease_assocs, temp_drugs, temp_druggability, temp_ab, temp_risk]:
                    remove(f)
                
                # Also, move/rename directories to something appropriate to the filename we're producing as a primary 
                # output
                html_reports_dir = sub("\.xlsx$", '_reports', self.outfile)
                gene_dir = html_reports_dir + "/" + gene_name
                if not isdir(html_reports_dir):
                    mkdir(html_reports_dir)
                if not isdir(gene_dir):
                    mkdir(gene_dir)
                move(report_html_file, gene_dir)
                move("figures", gene_dir)
                
                # self.progbar.emit((i / len(self.series)) * 100)
                # self.dialog_thread.update_progbar((i / len(self.series)) * 100)
                i += 1
        except Exception as e:
            print("Exception in report gen function:")
            print(e)
            print_exc()


class About(QDialog):
    def __init__(self):
        super(About, self).__init__()
        self.img = 'logo2.png'
        self.ui = about.Ui_About()
        self.ui.setupUi(self)
        self.setWindowIcon(QIcon('logo2.png'))
        self.ui.label_logo.setPixmap(QPixmap(getcwd() + '/' + self.img))


class BucketingHelp(QDialog):
    def __init__(self):
        super(BucketingHelp, self).__init__()
        self.ui = bucketing_help.Ui_BucketingHelp()
        self.ui.setupUi(self)
        self.setWindowIcon(QIcon('logo2.png'))
        invpad = 30
        self.ui.label_druggability_image_criteria.setPixmap(
            QPixmap(getcwd() + '/images/druggability.png').scaledToWidth(
                self.ui.label_druggability_image_criteria.width() - invpad))
        self.ui.label_druggability_image_tree.setPixmap(
            QPixmap(getcwd() + '/images/druggability_tree.png').scaledToWidth(
                self.ui.label_druggability_image_tree.width() - invpad))
        self.ui.label_safety_image_criteria.setPixmap(
            QPixmap(getcwd() + '/images/safety.png').scaledToWidth(
                self.ui.label_safety_image_criteria.width() - invpad))
        self.ui.label_safety_image_tree.setPixmap(
            QPixmap(getcwd() + '/images/safety_tree.png').scaledToWidth(
                self.ui.label_safety_image_tree.width() - invpad))
        self.ui.label_feasibility_image_criteria.setPixmap(
            QPixmap(getcwd() + '/images/feasibility.png').scaledToWidth(
                self.ui.label_feasibility_image_criteria.width() - invpad))
        self.ui.label_feasibility_image_tree.setPixmap(
            QPixmap(getcwd() + '/images/feasibility_tree.png').scaledToWidth(
                self.ui.label_feasibility_image_tree.width() - invpad))
        self.ui.label_antibodyability_image_criteria.setPixmap(
            QPixmap(getcwd() + '/images/antibodyability.png').scaledToWidth(
                self.ui.label_antibodyability_image_criteria.width() - invpad))
        self.ui.label_antibodyability_image_tree.setPixmap(
            QPixmap(getcwd() + '/images/antibodyability_tree.png').scaledToWidth(
                self.ui.label_antibodyability_image_tree.width() - invpad))
        self.ui.label_modality_image_criteria.setPixmap(
            QPixmap(getcwd() + '/images/modality.png').scaledToWidth(
                self.ui.label_modality_image_criteria.width() - invpad))
        self.ui.label_modality_image_tree.setPixmap(
            QPixmap(getcwd() + '/images/modality_tree.png').scaledToWidth(
                self.ui.label_modality_image_tree.width() - invpad))
        self.ui.label_tissue_engagement_image_criteria.setPixmap(
            QPixmap(getcwd() + '/images/tissue_engagement.png').scaledToWidth(
                self.ui.label_tissue_engagement_image_criteria.width() - invpad))
        self.ui.label_tissue_engagement_image_tree.setPixmap(
            QPixmap(getcwd() + '/images/tissue_engagement_tree.png').scaledToWidth(
                self.ui.label_tissue_engagement_image_tree.width() - invpad))


class Settings(QDialog):
    def __init__(self, api_keys):
        print(api_keys)
        super(Settings, self).__init__()
        self.ui = settings.Ui_SettingsDialog()
        self.ui.setupUi(self)
        self.setWindowIcon(QIcon('logo2.png'))
        if api_keys["OpenTargets_appName"]:
            self.ui.lineEdit_OpenTargets_appName.setText(api_keys["OpenTargets_appName"])
        if api_keys["OpenTargets_secret"]:
            self.ui.lineEdit_OpenTargets_secret.setText(api_keys["OpenTargets_secret"])
        self.ui.lineEdit_FDA_key.setText(api_keys["FDA_key"])
        self.ui.lineEdit_biogrid_key.setText(api_keys["biogrid_key"])
        self.ui.lineEdit_seq_similarity.setText(str(api_keys['seq_similarity_threshold']))
        self.ui.lineEdit_seq_similarity.setValidator(QIntValidator(1, 100, self))


class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.dir = getcwd()
        self.ui = gui.Ui_MainWindow()
        self.ui.setupUi(self)
        self.setWindowIcon(QIcon('logo2.png'))
        
        self.list_genes = None
        self.targetsfile = None
        self.infile = None
        self.outfile = None
        self.data = None
        self.list_of_targets = DataFrame()
        self.targetsmodel = None
        
        # Connect button pushed events to handler subs
        self.ui.pushButton_insertFromText.clicked.connect(self.insert_targets_from_text)
        self.ui.pushButton_insertFromFile.clicked.connect(self.insert_targets_from_file)
        self.ui.pushButton_clearList.clicked.connect(self.clear_list)
        self.data_downloaded = False
        self.change_download_button_role()
        self.ui.pushButton_download.setEnabled(False)
        self.ui.pushButton_openExistingData.clicked.connect(self.open_existing_data)
        self.ui.pushButton_set_disease_profile.clicked.connect(self.manage_disease_profiles)
        
        self.ui.actionScoring_guide.triggered.connect(self.launch_bucketing_help_dialog)
        self.ui.actionAbout.triggered.connect(self.launch_about_dialog)
        self.ui.actionImport_targets.triggered.connect(self.insert_targets_from_file)
        self.ui.actionOpen_existing_file.triggered.connect(self.open_existing_data)
        
        self.api_keys_path = "api_keys.json"
        self.api_keys = self.read_api_keys()
        self.ui.actionSettings.triggered.connect(self.launch_settings_dialog)
        
        self.downloaderThread = Downloader(api_keys=self.api_keys,
                                           literature=self.ui.checkBox_literature.isChecked(),
                                           ppi=self.ui.checkBox_ppi.isChecked())
        self.downloaderThread.progbar_update.connect(self.update_progbar)
        self.downloaderThread.status.connect(self.update_status)
        self.downloaderThread.got_data.connect(self.download_finished)
        self.downloaderThread.warnings.connect(self.warning_widget)
        self.ui.progressBar.setValue(0)
        
        self.current_disease_profile = None
        self.read_active_disease_profile()
        # Same for the 'about' dialog. 
        # self.ui.pushButton_about.clicked.connect(self.launch_about_dialog)
    
    def open_existing_data(self):
        self.update_status("Reading file...")
        self.infile = QFileDialog.getOpenFileName(self, 'Open existing druggability data', expanduser('~'),
                                                  "Excel files (*.xlsx)")[0]
        # Does file exist? If so, set self.outfile = None
        if isfile(self.infile):
            sleep(0.1)
            # Deactivate download/open-existing buttons!
            self.ui.pushButton_insertFromText.setEnabled(False)
            self.ui.pushButton_insertFromFile.setEnabled(False)
            self.outfile = None
            self.data_downloaded = True
            self.change_download_button_role()
            self.data = read_excel(self.infile, sheet_name=None)
            if 'External data input' in self.data.keys():
                # self.insert_targets_from_file(self.data['External data input'])
                self.update_input_targets_model(self.data['External data input'])
            self.launch_data_dialog()
        else:
            self.update_status("Ready")
            if self.infile:
                w = QWidget()
                QMessageBox.critical(w, "Error", "No data file found!\n" +
                                     "Selected file does not exist.'")
    
    def warning_widget(self, message):
        # Unbelievable that I have to try this, but here we are.
        w = QWidget()
        QMessageBox.critical(w, "Error", message)
    
    def update_progbar(self, i):
        self.ui.progressBar.setValue(i)
    
    def update_status(self, i):
        self.ui.lineEdit_status.setText(i)
    
    def read_active_disease_profile(self):
        path = "disease_profiles/active.txt"
        if exists(path):
            with open(path, "r") as f:
                name = f.read()
                if name:
                    print("Reading active disease profile " + name)
                    self.ui.lineEdit_disease_profile.setText(name)
                    profile = load(open("disease_profiles/" + name + ".json"))
                    self.current_disease_profile = profile
                    self.downloaderThread.disease_profile = profile
    
    def read_api_keys(self):
        api_keys = {"OpenTargets_appName": None,
                    "OpenTargets_secret": None,
                    "FDA_key": None,
                    "biogrid_key": None,
                    'seq_similarity_threshold': 40}
        if exists(self.api_keys_path):
            api_keys = load(open(self.api_keys_path))
        return api_keys
    
    def manage_disease_profiles(self):
        # self.update_status("Opening datasets required for disease profile specification...")
        # That doesn't work very well; better to pop up a dialog or something. 
        try:
            dlg = DiseaseProfileSelectorDialog(self.downloaderThread, self.current_disease_profile)
            if dlg.exec_():
                new_profile = dlg.make_some_json_output()
                self.ui.lineEdit_disease_profile.setText(new_profile['name'])
                # Save the name of this profile in disease_profiles/active.txt
                with open("disease_profiles/active.txt", "w") as f:
                    f.write(new_profile['name'])
                self.current_disease_profile = new_profile
                self.downloaderThread.disease_profile = new_profile
            self.update_status("Ready")
        except Exception as e:
            print("Exception in disease profile manager:")
            print(e)
            print_exc()
    
    def insert_targets_from_text(self):
        dlg = TextInsertDialog()
        if dlg.exec_():
            # If user input a list and clicked OK...
            new_targets, target_type = dlg.get_values()
            if len(new_targets) > 0:
                colnames = {'HGNC Symbol': 'HGNC Name',
                            'Ensembl ID': 'GeneID',
                            'Uniprot ID': 'Uniprot ID'}
                # Make a dataframe
                targets = DataFrame({colnames[target_type]: new_targets})
                # Add to existing/create new table model
                if self.targetsmodel:
                    targets = concat([self.targetsmodel.orig_data, targets], axis=0, sort=False)
                    targets = targets.drop_duplicates()
                targets['series'] = range(1, len(targets.index) + 1, 1)
                self.update_input_targets_model(targets)
    
    def update_input_targets_model(self, targets):
        # Put the 'series' column first
        targets = targets[['series'] + [i for i in targets.columns.tolist() if i != "series"]]
        if len(targets) > 0:
            self.ui.pushButton_download.setEnabled(True)
        self.targetsmodel = PandasModel(targets)
        self.ui.tableView.setModel(self.targetsmodel)
        self.ui.tableView.resizeColumnsToContents()
        self.list_of_targets = targets
        self.ui.label.setText(str(len(self.list_of_targets)) + " targets present")
        #                      Input a list of gene identifiers
    
    def insert_targets_from_file(self, existing_list=None):
        # We may be able to pass data here on opening an existing output file!
        # Act on it if so; otherwise, pick a file to open.
        if isinstance(existing_list, DataFrame):
            # Clear existing data, if any
            self.targetsmodel = None
            targets = existing_list
        else:
            self.targetsfile, filetype = QFileDialog.getOpenFileName(self,
                                                                     'Open existing druggability data', expanduser('~'),
                                                                     "Excel files (*.xlsx);;CSV files (*.csv);;All files (*)")
            
            # Good. Now, one way or another, we need to read that into pandas and have it show in our table view. 
            targets = None
            if filetype == "Excel files (*.xlsx)":
                targets = read_excel(self.targetsfile, index_col=None)
            elif filetype == "CSV files (*.csv)":
                targets = read_csv(self.targetsfile)
        
        if targets is not None:
            # Do some further validation here.
            # There are a few simple rules to apply. 
            # First, do any columns have 'GeneID' (any case) in? Make sure that's correctly cased.
            # May as well do series, HGNC name and Uniprot while I'm at it.
            for col in targets.columns.values.tolist():
                if col.upper() == 'GENEID' and col != 'GeneID':
                    targets = targets.rename(columns=lambda x: x.replace(col, 'GeneID'))
                if col.upper() in ["ENSEMBL ID", "ENSEMBL_ID", "ENSEMBL-ID", "ENSEMBLID"]:
                    targets = targets.rename(columns=lambda x: x.replace(col, 'GeneID'))
                if col.upper() == 'SERIES' and col != 'series':
                    targets = targets.rename(columns=lambda x: x.replace(col, 'series'))
                if col.upper() in ['HGNC NAME', 'HGNC_NAME', 'HGNC-NAME', 'HGNCNAME'] and col != 'HGNC Name':
                    targets = targets.rename(columns=lambda x: x.replace(col, 'HGNC Name'))
                if col.upper() in ['UNIPROT ID', 'UNIPROT_ID', 'UNIPROT-ID', 'UNIPROTID'] and col != 'Uniprot ID':
                    targets = targets.rename(columns=lambda x: x.replace(col, 'Uniprot ID'))
            if 'GeneID' not in targets.columns.values.tolist():
                # If not, do any of the column names have anything resembling 'ensembl ID' in them?
                ensembl_id_cols = [i for i in targets.columns.values.tolist()
                                   if search('ensembl[_\-\s.]id', i, IGNORECASE)]
                # If so, pick the first one with correctly formatted ensembl IDs and rename it to 'GeneID'
                for col in ensembl_id_cols:
                    test_value = targets[col].loc[targets[col].first_valid_index()]
                    if test_value:
                        # Is this an Ensembl ID?
                        if match('ENSG[0-9]{11}', test_value):
                            targets = targets.rename(columns=lambda x: x.replace(col, 'GeneID'))
                            break
            # Targets frame should now have at least one of GeneID, HGNC name, and Uniprot ID.
            required_cols = ['GeneID', 'HGNC Name', 'Uniprot ID']
            present_cols = list(set(targets.columns.values.tolist()).intersection(set(required_cols)))
            if present_cols:
                # There should also be values in at least one of those columns in every row. Drop any rows where that's
                # not the case - but do it non-destructively, for now.
                targets_na_dropped = targets.dropna(subset=present_cols, how="all")
                # If that did cause some rows to be dropped, ask the user what to do.
                diff = len(targets.index) - len(targets_na_dropped.index)
                if diff > 0:
                    msgBox = QMessageBox(self)
                    msgBox.setIcon(QMessageBox.Information)
                    msgBox.setText(str(diff) + " targets have no detectable identifiers.\n"
                                               "All targets must have at least one of an Ensembl ID,\n"
                                               "Uniprot ID or HGNC name.\n")
                    msgBox.setInformativeText("Proceed without them?")
                    msgBox.addButton(QMessageBox.Yes)
                    msgBox.addButton(QMessageBox.No)
                    
                    msgBox.setDefaultButton(QMessageBox.No)
                    ret = msgBox.exec_()
                    
                    if ret == QMessageBox.Yes:
                        targets = targets_na_dropped
                        print("Dropped " + str(diff) + " targets with no identifiers")
                    else:
                        return
                self.data_downloaded = False
                self.change_download_button_role()
                # Get existing data, if any, to merge with new data
                if self.targetsmodel:
                    targets = concat([self.targetsmodel.orig_data, targets], axis=0, sort=False)
                    targets = targets.drop_duplicates()
                
                # Overwrite the series column (whether or not it exists); this is super important for several things.
                targets['series'] = range(1, len(targets.index) + 1, 1)
                
                self.update_input_targets_model(targets)
            else:
                # If columns are missing, give a warning message to that effect.
                w = QWidget()
                QMessageBox.critical(w, "Error", "Insufficient data found!\n" +
                                     "Targets must have at least one of an Ensembl ID,\n"
                                     "Uniprot ID or HGNC name.'")
    
    def clear_list(self):
        # These two buttons will have been disabled on clicking the 'Download' button. Enable them here.
        if not self.ui.pushButton_insertFromText.isEnabled():
            self.ui.pushButton_insertFromText.setEnabled(True)
        if not self.ui.pushButton_insertFromFile.isEnabled():
            self.ui.pushButton_insertFromFile.setEnabled(True)
        if self.targetsmodel:
            self.targetsmodel.delete_everything()
        if self.ui.pushButton_download.text() == "View downloaded\ndruggability data":
            self.data_downloaded = False
            self.change_download_button_role()
        self.ui.pushButton_download.setEnabled(False)
        self.ui.label.setText("Input a list of gene identifiers")
    
    def download_data(self):
        # set prog bar to 0%.
        self.ui.progressBar.setValue(0)
        # Pick an output file
        default_filename = expanduser('~') + "\\" + "druggability.xlsx"
        if self.targetsfile:
            default_filename = expanduser('~') + "\\" + sub('.xlsx|.csv', '',
                                                            basename(self.targetsfile)) + "_druggability.xlsx"
        self.outfile = QFileDialog.getSaveFileName(self, 'Save output as...', default_filename,
                                                   "Excel files (*.xlsx)")[0]
        # Now actually fire it up and start downloadin' stuff
        # (Only do that if an output file has been set!)
        if self.outfile:
            self.downloaderThread.set_target_list(self.list_of_targets)
            self.downloaderThread.start()
            # Deactivate download/open-existing buttons!
            self.ui.pushButton_insertFromText.setEnabled(False)
            self.ui.pushButton_insertFromFile.setEnabled(False)
            self.ui.pushButton_clearList.setEnabled(False)
            self.ui.pushButton_set_disease_profile.setEnabled(False)
            self.ui.pushButton_download.setEnabled(False)
            self.ui.pushButton_openExistingData.setEnabled(False)
            self.ui.checkBox_ppi.setEnabled(False)
            self.ui.checkBox_literature.setEnabled(False)
    
    def download_finished(self, data):
        self.data = data
        self.data_downloaded = True
        self.change_download_button_role()
        self.update_status("Saving downloaded data...")
        # Write it to the file, if one is specified
        self.write_data(self.outfile, self.data)
        # Reactivate these guys now that download is complete
        # self.ui.pushButton_insertFromText.setEnabled(True)
        # self.ui.pushButton_insertFromFile.setEnabled(True)
        self.ui.pushButton_clearList.setEnabled(True)
        self.ui.pushButton_set_disease_profile.setEnabled(True)
        self.ui.pushButton_download.setEnabled(True)
        self.ui.pushButton_openExistingData.setEnabled(True)
        self.ui.checkBox_ppi.setEnabled(True)
        self.ui.checkBox_literature.setEnabled(True)
        # Need to send this data off to the window that displays it in a nice sortable, searchable tab format.
        self.update_status("Opening data viewer...")
        self.launch_data_dialog()
    
    def write_data(self, filename, sheets):
        try:
            writer = ExcelWriter(filename)
            for name in sheets:
                df = sheets[name]
                df.to_excel(writer, name, index=False)
                # This cunning approach mostly comes from https://stackoverflow.com/a/40535454
                worksheet = writer.sheets[name]  # pull worksheet object
                for idx, col in enumerate(df):  # loop through all columns
                    series = df[col]
                    max_len = max((
                        series.astype(str).map(len).max(),  # len of largest item
                        len(str(series.name))
                        # len of column name/header (which may be substantially longer than any data values)
                    )) + 1  # adding a little extra space
                    # Put a hard cap on that, or we may get huge, unmanageably wide columns.
                    if max_len > 85:
                        max_len = 85
                    worksheet.set_column(idx, idx, max_len)  # set column width
            writer.save()
        except Exception as e:
            w = QWidget()
            QMessageBox.critical(w, "Error", str(e) + "\nClick OK to retry.")
            self.write_data(filename, sheets)
    
    def launch_data_dialog(self):
        # data can be either a list of dataframes or a path to an excel file
        print("Launching data window...")
        datawin = DataWindow(self.data)
        datawin.make_reports.connect(self.launch_report_generator)
        datawin.show()
        datawin.exec_()
        self.update_status("Ready")
    
    def launch_report_generator(self, make_report_for_these):
        print("Launching report generator thing...")
        sourcefile = self.outfile if self.outfile else self.infile
        # It turns out this HAS to be run from the main thread, in order to work at all. 
        repwin = ReportGenerator(df=self.data,
                                 series=make_report_for_these,
                                 outfile=sourcefile)
        repwin.make_reports()
        print("Reports should be done. Check the directory where output got saved.")
    
    def change_download_button_role(self):
        self.ui.pushButton_download.disconnect()
        if self.data_downloaded:
            self.ui.pushButton_download.setText("View downloaded\ndruggability data")
            self.ui.pushButton_download.clicked.connect(self.launch_data_dialog)
        else:
            self.ui.pushButton_download.setText("Download\ndruggability data")
            self.ui.pushButton_download.clicked.connect(self.download_data)
    
    def launch_about_dialog(self):
        about = About()
        about.show()
        about.exec_()
    
    def launch_bucketing_help_dialog(self):
        b_help = BucketingHelp()
        b_help.show()
        b_help.exec_()
    
    def launch_settings_dialog(self):
        sett = Settings(self.api_keys)
        if sett.exec_():
            self.api_keys = {"OpenTargets_appName": sett.ui.lineEdit_OpenTargets_appName.text(),
                             "OpenTargets_secret": sett.ui.lineEdit_OpenTargets_secret.text(),
                             "FDA_key": sett.ui.lineEdit_FDA_key.text(),
                             "biogrid_key": sett.ui.lineEdit_biogrid_key.text(),
                             "seq_similarity_threshold": int(sett.ui.lineEdit_seq_similarity.text())}
            self.downloaderThread.set_api_keys(self.api_keys)
            with open(self.api_keys_path, 'w') as outfile:
                dump(self.api_keys, outfile)


if __name__ == "__main__":
    app = QApplication(argv)
    
    # Splash screen
    splash_pix = QPixmap('splash.jpg')
    splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    splash.show()
    sleep(1)
    
    mw = MainWindow()
    mw.show()
    splash.finish(mw)
    exit(app.exec_())
