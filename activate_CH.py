import sys
import pickle
from rdkit import Chem
from rdkit.Chem import AllChem

filename = sys.argv[1]
rxn = sys.argv[2]

rm_proton = ( ('[CX4;H3]','[CH2-]'),('[CX4;H2]','[CH-]'),('[CX4;H1]','[C-]'),('[CX3;H2]','[CH-]'),
              ('[CX3;H1]','[C-]'), ('[CX2;H1]','[C-]'))

rm_hydride =( ('[CX4;H3]','[CH2+]'),('[CX4;H2]','[CH+]'),('[CX4;H1]','[C+]'),('[CX3;H2]','[CH+]'),
              ('[CX3;H1]','[C+]'), ('[CX2;H1]','[C+]'))

rm_hydrogen = ( ('[CX4;H3]','[CH2]'),('[CX4;H2]','[CH]'),('[CX4;H1]','[C]'),('[CX3;H2]','[CH]'),
              ('[CX3;H1]','[C]'), ('[CX2;H1]','[C]'))

if rxn == "rm_proton":
     smartsref = rm_proton
     delta_charge = -1
if rxn == "rm_hydride":
     smartsref = rm_hydride
     delta_charge = 1
if rxn == "rm_hydrogen":
     smartsref = rm_hydrogen
     delta_charge = 0

def check_chiral(m,ion_smiles):
    chiral_atoms_mol = Chem.FindMolChiralCenters(m)
    number_of_chiral_atoms_mol = len(chiral_atoms_mol)
    if number_of_chiral_atoms_mol == 0:
        return ion_smiles

    ion = Chem.MolFromSmiles(ion_smiles)
    chiral_atoms_ion = Chem.FindMolChiralCenters(ion)

    for chiral_atom_mol,chiral_atom_ion in zip(chiral_atoms_mol,chiral_atoms_ion):
        if chiral_atom_mol != chiral_atom_ion:
            ion.GetAtomWithIdx(chiral_atom_mol[0]).InvertChirality()
        ion_smiles = Chem.MolToSmiles(ion,isomericSmiles=True)
#       print Chem.FindMolChiralCenters(ion)

    return ion_smiles

def change_mol(name, charge, mol, smartsref):
    for (smarts1, smiles2) in smartsref:
        patt1 = Chem.MolFromSmarts(smarts1)
        patt2 = Chem.MolFromSmiles(smiles2)
        if(mol.HasSubstructMatch(patt1)):
            atoms = mol.GetSubstructMatches(patt1)
#convert tuple of tuple to list: ((1,),(2,)) -> [1,2]
            atoms = [element for tupl in atoms for element in tupl]
            newmol = AllChem.ReplaceSubstructs(mol, patt1, patt2)
            for atom,ion in zip(atoms,newmol):
                ion_smiles = Chem.MolToSmiles(ion,isomericSmiles=True)
                ion_smiles = check_chiral(mol,ion_smiles)
                if delta_charge != 0:
                    new_charge = str(charge+delta_charge)
                else:
                    new_charge = "rad"

                print name+"_"+new_charge+"="+str(atom), ion_smiles

    return


with open(filename, "r") as file:

    for line in file:
        name,smiles = line.split()
        mol = Chem.MolFromSmiles(smiles)
        Chem.Kekulize(mol,clearAromaticFlags=True)
        charge = Chem.GetFormalCharge(mol)
        print name+"_"+str(charge)+"=0", smiles
    
        change_mol(name, charge, mol, smartsref)

