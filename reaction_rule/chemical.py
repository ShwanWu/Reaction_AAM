import os
import glob
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdChemReactions import ChemicalReaction
from rdkit.Chem import Draw

class Chemical:

    def __init__(self, mol):
        # self is a mol object
        self.mol = mol

    def get_atom_feature(self, atom):
        """
        get one atom feature
        :param: the atom and the molecule that the atom belongs
        :return: a dict represent feature of one atom, use neighbours'
                mapped number(identity because every atom has different number)) and bond type
                {'neigh_map': [the mapped_number of neighbors],
                 'neigh_bond': [the bond type(int) of neighbors]}
        """
        idx = atom.GetIdx()
        feat = dict()
        # old_list is the list that before sorted.
        old_map_list = []
        old_bond_list = []
        for neigh in atom.GetNeighbors():
            idn = neigh.GetIdx()
            neigh_map_num = neigh.GetAtomMapNum()
            neigh_bond_type = int(self.GetBondBetweenAtoms(idx, idn).GetBondType())
            old_map_list.append(neigh_map_num)
            old_bond_list.append(neigh_bond_type)
        # return a list of index of components in map_list
        # the index list is sort by number
        # eg: mapped_list = [12, 5, 7] => index_list = [1, 2, 0]
        # where 1 means the index 1 (refer 5 in mapped_list) should be placed first order
        index = np.argsort(old_map_list)
        length = len(old_bond_list)
        # new_list is the old list after sorted
        new_map_list = []
        new_bond_list = []
        for i in range(length):
            order = index[i]
            new_map_list.append(old_map_list[order])
            new_bond_list.append(old_bond_list[order])
        feat['neigh_map'] = new_map_list
        feat['neigh_bond'] = new_bond_list
        return feat

    def get_mol_feature(self):
        """
        get feature collection of a molecule
        :param: mol object of a molecule
        :return: feature_list of a molecule contain every atom's feature
        """
        neigh_map = dict()  # it is used to extend changed atoms to radius 1
        mol_feat = []
        for atom in self.GetAtoms():
            AAM = atom.GetAtomMapNum()
            atom_feat = Chemical.get_atom_feature(self, atom)
            mol_feat.append([AAM, atom_feat])
            neigh_map[AAM] = atom_feat['neigh_map']
        return mol_feat, neigh_map

    def find_changed_atoms(self, other_mol):
        """
        find the feature changed atom of self.mol compare with another_mol
        :param: self.mol other_mol
        :return: a list contain list of changed atoms
        """
        feat = Chemical.get_mol_feature(self)[0]
        other_feat = Chemical.get_mol_feature(other_mol)[0]
        changed_atoms = []
        for a in feat:
            # initial exit as False because we haven't check all the atoms
            exit = False
            for b in other_feat:
                # if the same atoms AAM are found, then compare their feature
                if a[0] == b[0]:
                    exit = True
                    if a[1] != b[1]:
                        changed_atoms.append(a[0])
            if exit == False:
                changed_atoms.append(a[0])
        return changed_atoms

    def extend_react_atoms_r(self, atom_list, r):
        pass

    def show_image(self):
        """
        show the image representation of molecule of a reaction
        :param: a list of mol object of one reaction
        """
        img = Draw.MolsToGridImage([self])
        img.show()
