# Convert chemist-friendly SDF to programmer-friendly CSV.
from rdkit.Chem import PandasTools

SDFFile = "Data/withdrawn_withdrawn_compounds.sdf"
BRDLigs = PandasTools.LoadSDF(SDFFile)
BRDLigs.to_csv("Data/withdrawn_withdrawn_compounds.csv", index=False)
