import sys
import math
import subprocess
sys.path.append("/usr/local/lib/python2.7/site-packages/")
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
#import small_ref_set as rf


filename = sys.argv[1]
filename2 = sys.argv[2]
filename3 = "reference_"+filename.split(".")[0].split("_",1)[1]+".energies"
if len(sys.argv) > 2:
   pKa_cutoff = float(sys.argv[3])
else:
   pKa_cutoff = 1.0


delta_charge = 0
if "proton" in filename:
    delta_charge = -1

if "hydride" in filename:
    delta_charge = 1

if "hydrogen" in filename:
    delta_charge = 0

energies = {}
log_files = {}
names = []
ions = []
ion2smiles = {}
ref_energies = {}
atom_radius = 2

smarts_ref = ( ('[C;H3]',10.4,'aliphatic'),
               ('[C;H2]',10.4,'aliphatic'),
               ('[C;H1]',10.4,'aliphatic'),
               ('[c;H1]',10.6,'aromatic') )

class pKa_class:
    def __init__(self,pKa,deprotonated_atom,prot_file):
        self.pKa = pKa
        self.deprotonated_atom = deprotonated_atom
        self.prot_file = prot_file

def find_ref(reactant_name,product_name):

    deprotonated_atom = product_name.split("=")[1]
    deprotonated_atom = int(deprotonated_atom)
    reactant_smiles = ion2smiles[reactant_name] 
    product_smiles = ion2smiles[product_name] 

    reactant = Chem.MolFromSmiles(reactant_smiles)
    product = Chem.MolFromSmiles(product_smiles)

    max_length = 0
    nref = "none"
    pKa_nref = 0.0
    for smarts, pKa, name in smarts_ref:
        ref_mol = Chem.MolFromSmarts(smarts)

        matches_tuple = reactant.GetSubstructMatches(ref_mol) 
        matches_list = [element for tupl in matches_tuple for element in tupl]
        if reactant.HasSubstructMatch(ref_mol) and deprotonated_atom in matches_list:
           length = len(list(reactant.GetSubstructMatches(ref_mol)))
           if length > max_length:
              nref = name
              max_length = length
              pKa_nref = pKa
    
    if nref == "aliphatic":
       nref = "toluene_0"

    if nref == "aromatic":
       nref = "toluene_0"

    return nref, deprotonated_atom, pKa_nref, reactant_smiles, product_smiles;

def find_better_ref(reactant,deprot_atom, ref_list):
    nref = "none"
    pKa_nref = 0.0
    max_length = 0
    for smiles, pKa, name in ref_list:
        ref_mol = Chem.MolFromSmiles(smiles)
#       ref_mol = Chem.MolFromSmarts(smiles)
 
    #    print Chem.MolToSmiles(reactant), smiles,deprot_atom, list(reactant.GetSubstructMatch(ref_mol))

        if reactant.HasSubstructMatch(ref_mol) and deprot_atom in list(reactant.GetSubstructMatch(ref_mol)):
           length = len(list(reactant.GetSubstructMatch(ref_mol)))
           if length > max_length:
              nref = name
              max_length = length
              pKa_nref = pKa

    return nref, pKa_nref

with open(filename3, "r") as reference_energies_file:
    for line in reference_energies_file:
        words = line.split()
        log_file  = words[0]
        ion = words[0].split("=")[0]
        if words[1] == "CURRENT":
            energy = float(words[-1])
        else:
            energy = float(words[6])
        if ion not in ref_energies:
           ref_energies[ion] = energy
        elif energy < ref_energies[ion]:
           ref_energies[ion] = energy

with open(filename2, "r") as smiles_file:
    for line in smiles_file:
        words = line.split()
        ion = words[0]
        smiles = words[1]
        ion2smiles[ion] = smiles 
    
with open(filename, "r") as energies_file:
    for line in energies_file:
        words = line.split()
        log_file  = words[0]
        ion = words[0].split("+")[0]
        name = ion.split("=")[0]
        if words[1] == "CURRENT":
            energy = float(words[-1])
        else:
            energy = float(words[6])
        if name not in names:
            names.append(name)
        if ion not in ions:
            ions.append(ion)
        if ion not in energies:
           energies[ion] = energy
           log_files[ion] = log_file
        elif energy < energies[ion]:
           energies[ion] = energy
           log_files[ion] = log_file

for name in names:
    if delta_charge != 0:
        charge = int(name.split("_")[1])    
        name_of_deprotonated = name.split("_")[0]+"_"+str(charge+delta_charge)
    else:
        if name.split("_")[1] == "0":
           name_of_deprotonated = name.split("_")[0]+"_rad"
        else:
           name_of_deprotonated = "skip"

    if name_of_deprotonated not in names:
        continue

# find ions containing name and name_of_deprotonated and collect in lists
    protonated = [item for item in ions if name in item]
    deprotonated = [item for item in ions if name_of_deprotonated in item]

    pKa_objects = []
    for deprot in deprotonated:
        for prot in protonated:
            ref,deprotonated_atom, pKa_ref, r_smiles, p_smiles = find_ref(prot,deprot)
            if ref == "skip": 
                continue
            delta_G = energies[prot] - energies[deprot]
            if delta_charge != 0:
                ref_deprot = ref.split("_")[0]+"_"+str(int(ref.split("_")[1])+delta_charge)
            else:
                ref_deprot = ref.split("_")[0]+"_rad"
            delta_G_ref = ref_energies[ref] - ref_energies[ref_deprot]
            pKa = pKa_ref + (delta_G_ref - delta_G)/1.36
#           print prot, deprot, deprotonated_atom, pKa
            pKa_object = pKa_class(pKa,str(deprotonated_atom),log_files[deprot])
            pKa_objects.append(pKa_object)

# sort objects from lowest to highest pKa value
    pKa_objects.sort(key=lambda x: x.pKa)
    pKa_min = pKa_objects[0].pKa
    deprot_atoms = []
    prot_files = []
# collect atoms with pKa values within cutoff of lowest pKa value
    for object in pKa_objects:
        if object.pKa - pKa_min <= pKa_cutoff:
            deprot_atoms.append(object.deprotonated_atom)
            prot_files.append(object.prot_file)

    print name.split("_")[0], ",".join(deprot_atoms), ",".join(prot_files)
