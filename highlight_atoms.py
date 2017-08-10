import os,sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole

name2smiles = {}

pka_file_name = sys.argv[1]
smiles_file_name = pka_file_name.split("_")[0]+".smiles"
svg_file_name = pka_file_name.split(".")[0]+".svg"



with open(smiles_file_name, "r") as smiles_file:
    for line in smiles_file:
        name, smiles = line.split()
        name2smiles[name] = smiles

mols = []
names = []
atoms = []

with open(pka_file_name, "r") as pka_file:
    for line in pka_file:
        words = line.split()
        name = words[0]
        smiles = name2smiles[name]
        atoms.append(map(int,words[1].split(",")))
        m = Chem.MolFromSmiles(smiles)
        mols.append(m)
        Chem.Kekulize(m)
        names.append(name)
        Draw.DrawingOptions.includeAtomNumbers=True

img = Draw.MolsToGridImage(mols,molsPerRow=4,subImgSize=(200,200),legends=[x for x in names],useSVG=True,highlightAtomLists=atoms)


print svg_file_name
with open(svg_file_name, 'w') as svg_file:
    svg_file.write(img.data)

os.system('sed -i "s/xmlns:svg/xmlns/" '+svg_file_name)
