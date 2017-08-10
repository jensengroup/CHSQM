import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

filename = sys.argv[1]

file = open(filename, "r")

for line in file:
    words = line.split()
    name = words[0]
    smiles = words[1]

    m = Chem.AddHs(Chem.MolFromSmiles(smiles))

    rot_bond = rdMolDescriptors.CalcNumRotatableBonds(m)
    print name, rot_bond
    confs = min(1 + 3*rot_bond,20)
#   rot_bond = rdMolDescriptors.CalcNumRotatableBonds(Chem.MolFromSmiles(smiles))
#   print name, rot_bond

    AllChem.EmbedMultipleConfs(m,numConfs=confs,useExpTorsionAnglePrefs=True,useBasicKnowledge=True)
#   AllChem.EmbedMultipleConfs(m,numConfs=confs)
#   AllChem.MMFFOptimizeMoleculeConfs(m,maxIters=1000)

    energies = []
    for i,conf in enumerate(m.GetConformers()):
        tm = Chem.Mol(m,False,conf.GetId())
        file = name+"+"+str(i)+".sdf"
        writer = Chem.SDWriter(file)
        writer.write(tm)
