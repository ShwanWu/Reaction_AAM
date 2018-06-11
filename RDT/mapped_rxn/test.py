from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdChemReactions
from rdkit.Chem.rdChemReactions import ChemicalReaction
import os

# reaction_smiles = "CC(O)CC(=O)OC(C)CC(O)=O.O>>OC(=O)CC(C)O.CC(O)CC(O)=O"
# cmd1 = "cd Documents/"
# os.system(cmd1)
# cmd2 = "java -jar ReactionDecoder.jar -Q SMI -q " + "reaction_smiles" + " -g -j AAM -f TEXT"
# os.system(cmd2)


file_path = '/home/ywu672/Documents/bitbucket/RDT/mapped_rxn/ECBLAST_smiles_AAM.rxn'
rxn = AllChem.ReactionFromRxnFile(file_path)
print(rxn)
ChemicalReaction.Initialize(rxn)
print(ChemicalReaction.GetReactants(rxn))
changed_atoms = ChemicalReaction.GetReactingAtoms(rxn)
print(changed_atoms)

mol = ChemicalReaction.GetReactantTemplate(rxn, 0)
print(mol)

from rdkit.Chem import Draw

img = Draw.MolsToGridImage([mol])
img.show()
