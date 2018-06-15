import os
import glob
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Chem import Draw
from chemical import Chemical

'''
use RDT to read a reaction smiles and generate rxn file, 
then find the reaction rule with user_defined radius form center changed atom.

:param  SMILES_input: a string of reaction SMILES read from database
:param  code_path: the path where this code is placed
:param  radius: the radius that the reaction core is limited to
'''

def delete_his(path):
    for infile in glob.glob(os.path.join(path, '*.log')):
        os.remove(infile)

def get_rxn(file_path):
    """
    get ChemicalReaction object from file and initialize
    :param: file_path of rxn
    :return: rxn object
    """
    rxn = AllChem.ReactionFromRxnFile(file_path)
    ChemicalReaction.Initialize(rxn)
    return rxn

def get_reactant(rxn, int):
    """
    get reactant mol object
    :param: rxn object and reactant index
    :return: mol object of reactant
    """
    mol = ChemicalReaction.GetReactantTemplate(rxn, int)
    return mol

def get_product(rxn, int):
    """
    get product mol object
    :param: rxn object and product index
    :return: mol object of product
    """
    mol = ChemicalReaction.GetProductTemplate(rxn, int)
    return mol

def get_AAM_dict(mol):
    """
    For a molecule, the idx assign b is different with the AAM
    In order to get bond between two know AAM atoms, we need use idx as connection
    :param mol: molecule mol object
    :return: a dict:{key = AAM, value = idx}
    """
    AAM_dict = dict()
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        AAM = atom.GetAtomMapNum()
        AAM_dict[AAM] = idx
    return AAM_dict

def get_map_list(mol):
    """get a list of mapped number of a molecule,
    basic used for find if one number is in the this molecule or not
    :param  molecule mol object
    :return a list = [AAM of each atom]
    """
    mapped_list = []
    for atom in mol.GetAtoms():
        mapped_list.append(atom.GetAtomMapNum())
    return mapped_list



def find_react_atoms(reactant_list, product):
    """
    find the feature changed atom of reactant of product
    :param: reactants, product in mol object
    :return: a list contain list of changed atoms in reactant and list of changed atoms in product
            [[[changed atom in R1],[changed atom in R2]], [changed atoms in P1]]
    """
    # create list storing the changed atoms of each reactant
    R_atom_list = []
    P_atom_list = []
    for i in range(len(reactant_list)):
        r_atom_list = Chemical.find_changed_atoms(reactant_list[i], product)
        R_atom_list.append(r_atom_list)
        p_atom_list = Chemical.find_changed_atoms(product, reactant_list[i])
        P_atom_list.append(p_atom_list)
    # for product changed atoms, we need find the union of every changed list when compare with every reactant
    new_P_atom_list = []
    for j in range(len(P_atom_list)):
        new_P_atom_list = list(set(new_P_atom_list).union(set(P_atom_list[j])))
    return R_atom_list, new_P_atom_list

def extend_react_atoms_e(reactant_list, product):
    """
    extent the changed atom by environment rules
    :param: list of changed atoms of reactants and product
    :return: reaction core atoms list
    """
    R_center, P_center = find_react_atoms(reactant_list, product)  # [[4, 5], [19]], [4, 5, 19]
    # initialize a copy of center atoms.
    new_R_center = R_center[:]
    new_P_center = P_center[:]

    for i in range(len()):
        pass
    return new_R_center, new_P_center

def extend_react_atoms_r(reactant_list, product, r):
    """
    extent the changed atom to neigh with radius 1
    :param: list of changed atoms of reactants and product
    :param: r
    :return: same form list after extended
    """
    R_center, P_center = find_react_atoms(reactant_list, product)  # [[4, 5], [19]], [4, 5, 19]
    # initialize a copy of center atoms.
    new_R_center = R_center[:]
    new_P_center = P_center[:]
    # use n to control iteration, every 'while' use former center atoms as basic atoms to extend
    n = 0
    while n < r:
        n += 1
        # every time when union the old center and neighbour atoms, update the center atom.
        for i in range(len(R_center)):
            # for every reactant changed atom list
            # get a dict: atom AAM with neighbor AAM. eg: {1:[2], 2:[1,3,14], 14:[2], 3:[2,4]}
            neigh_map = Chemical.get_mol_feature(reactant_list[i])[1]
            for j in range(len(R_center[i])):
                # for every mapped number, find the neighbor(radius=1) atoms' mapped number
                # we can find the neighbor atoms' mapped number from it's atom_feature
                center_atom_AAM = R_center[i][j]
                neighbor_atom_AAM = neigh_map[center_atom_AAM]
                new_R_center[i] = list(set(new_R_center[i]).union(set(neighbor_atom_AAM)))
        R_center = new_R_center[:]

        for k in range(len(P_center)):
            center_atom_AAM = P_center[k]
            neighbor_atom_AAM = Chemical.get_mol_feature(product)[1][center_atom_AAM]
            new_P_center = list(set(new_P_center).union(set(neighbor_atom_AAM)))
        P_center = new_P_center[:]

    return new_R_center, new_P_center

def get_reaction_core(reactant_list, product, n):
    """
    get the reaction_core in a mol representation, which is consist of extended changed atoms
    first find the bonds that in the reaction core  as connection and then construct a new mol
    :param: all reactants and products
    :return a mol object of the reaction core(substructure)
    """
    R_center, P_center = extend_react_atoms_r(reactant_list, product, n)
    R_core_smiles = []
    P_core_smiles = []
    # for every reactants:
    for i in range(len(reactant_list)):
        # find all the bonds between the atoms in the list like AAM=[1,2,3,14] => idx=[x,x,x,x]
        list = R_center[i]
        dict = get_AAM_dict(reactant_list[i])
        atoms_to_use = []
        # limit the atom symbol to generate SMILES
        atom_symbol_list = [atom.GetSymbol() for atom in reactant_list[i].GetAtoms()]
        for j in list:
            idx = dict[j]
            atoms_to_use.append(idx)
        smiles = Chem.MolFragmentToSmiles(reactant_list[i], atoms_to_use, atomSymbols=atom_symbol_list)
        R_core_smiles.append(smiles)
    list = P_center
    atoms_to_use = []
    dict = get_AAM_dict(product)
    atom_symbol_list = [atom.GetSymbol() for atom in product.GetAtoms()]
    for j in list:
        idx = dict[j]
        atoms_to_use.append(idx)
    smiles = Chem.MolFragmentToSmiles(product, atoms_to_use, atomSymbols=atom_symbol_list)
    P_core_smiles.append(smiles)
    return R_core_smiles, P_core_smiles

def get_reaction_rule(reactant_list, product, n):
    """
    use lists of fragments to generate reaction rule
    :param reactant_list: a list of reactants core in smiles
    :param product: a list of product core in smiles
    :return: final reaction rule in SMILES
    """
    R_core_smiles, P_core_smiles = get_reaction_core(reactant_list, product, n)
    SMILES = ""
    for i in range(len(R_core_smiles)):
        if i == 0:
            SMILES += str(R_core_smiles[i])
        else:
            SMILES += '.' + str(R_core_smiles[i])
    SMILES += '>>'
    for j in range(len(P_core_smiles)):
        if j == 0:
            SMILES += str(P_core_smiles[j])
        else:
            SMILES += '.' + str(P_core_smiles[j])
    return SMILES

if __name__ == '__main__':


    code_path = "/home/ywu672/Documents/bitbucket/Reaction_rule/AAM/reaction_rule"

    SMILES_input = "FC1=CC=C(C=C[N+](=O)[O-])C=C1.C1(CCCCC1)=O" \
             ">>FC1=CC=C(C=C1)[C@H](C[N+](=O)[O-])[C@H]1C(CCCC1)=O"

    radius = 1

    input = "\"" + SMILES_input + "\""
    print(SMILES_input)
    cmd = 'java -jar ReactionDecoder.jar -Q SMI -q ' + input + ' -g -j AAM -f TEXT'
    os.system(cmd)
    delete_his(code_path)

    file_path = code_path + '/ECBLAST_smiles_AAM.rxn'
    rxn = get_rxn(file_path)

    # get rdkit.Chem.rdchem.Mol object of Reactants and products
    R1 = get_reactant(rxn, 0)
    R2 = get_reactant(rxn, 1)
    P1 = get_product(rxn, 0)
    R_list  = [R1, R2]
    mols = [R1, R2, P1]

    #find changed atoms
    print(find_react_atoms(R_list, P1))
    # #find changed atoms
    # print(extend_react_atoms_r(R_list, P1, radius))
    # # find list of reaction core
    # print(get_reaction_core(R_list, P1, radius))
    # # find reaction rule
    # print(get_reaction_rule(R_list, P1, radius))


