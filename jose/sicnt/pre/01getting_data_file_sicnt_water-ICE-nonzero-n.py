import collections 
from collections import Counter
import pandas as pd
import numpy as np

def main():
    nanchors = 3
    density = 0.471
    charge_si = 0.6
    atoms_per_anchor = 24
    with open('final_pdb_{}.pdb'.format(density), 'r') as infile, open('data-file-sicnt-water-{}.data'.format(charge_si), 'w') as outfile, \
        open('silicon.txt', 'w') as silicon_out,  open('carbon.txt', 'w') as carbon_out, open('hydrogen.txt', 'w') as hydro_out, \
        open('oxygen.txt', 'w') as oxy_out, open('water.txt', 'w') as water_out,  open('sicnt_no_anchors.txt', 'w') as sicar_no_anch_out, \
        open('sicnt.txt', 'w') as sicnt_out, open('anchors_T.txt', 'w') as anchors_sicnt_T, \
        open('almost_all_ids.txt', 'w') as sicnt_no_anchors_and_water_out, \
        open('sicnt_gr_1.txt', 'w') as sicnt_group1, open('sicnt_gr_2.txt', 'w') as sicnt_group2, \
        open('nofrozen_molecules.txt', 'r') as nofroz, open('frozen_molecules.txt', 'r') as froz:
        
        nofrozm = []
        frozm = []
        
        for i in nofroz:
            il = i.strip()
            nofrozm.append(il)
        for j in froz:
            frozm.append(j)
            
        print(nofrozm)
        
        
        #masses: 
        #Si(1) = 28.0850
        #C(2) = 12.0110
        #H(3) =  1.0080
        #O(4)  = 15.9994

        atom_list, conections, silicon, carbon, hydrogen, oxygen, sicnt  = ([] for i in range(7))
        num_atoms = 0
        for j, lines in enumerate(infile):
            if j == 2:
                di = lines.split()
                dimen = [float(di[1]), float(di[2]), float(di[3])]    
            
            elif j > 8:    
                elementspdb = lines.split()
                
                if elementspdb[0] =='ATOM':
                    
                    #get atomic positions. "Atoms" part of data file
                    atomic_list = atom_line(elementspdb[1], elementspdb[5], elementspdb[6], elementspdb[7], elementspdb[10])
                    atom_list.append(atomic_list.get_q_type_mol_tag())
                    
                    #get groups of atoms, oxygen, carbon, hydrogen and water to be used in input file
                    atom_types = split_types(elementspdb[1], elementspdb[10])
                    
                    if atom_types.get_Si() != None:
                        silicon.append(atom_types.get_Si())
                    
                    elif atom_types.get_C() != None:
                        carbon.append(atom_types.get_C())
                        
                    elif atom_types.get_H() != None:
                        hydrogen.append(atom_types.get_H())
                        
                    elif atom_types.get_O() != None:
                        oxygen.append(atom_types.get_O())
                        
                    num_atoms += 1
                        
                elif elementspdb[0] =='CONECT':
                    conections.append(elementspdb[1:])
        
        
        print(*silicon, file = silicon_out)
        print(*carbon, file = carbon_out)
        print(*hydrogen, file = hydro_out)
        print(*oxygen, file = oxy_out)
        water = hydrogen + oxygen
        water_s = sorted(water)
        sicnt = silicon + carbon
        sicnt_s = sorted(sicnt)
        print(*water_s, file = water_out)
        print(*sicnt_s, file = sicnt_out)
        
        #for l in oxygen:
            #print(l)
        
        
        atom_list_obj = adding_numeration_to_items(atom_list)#
        atom_types_dic = atom_list_obj.count_atom_types()
        #print(atom_types_dic)
    
        all_bonds = getting_bonds(conections, carbon, silicon, hydrogen, oxygen)
        all_angles = getting_angles(conections, carbon, hydrogen, oxygen)

        water_items = getting_water_items(all_bonds, all_angles)
        water_bonds = water_items.get_water_bonds()
        water_angles = water_items.get_water_angles()
        print(all_bonds)
        print(water_bonds)
        #print(all_angles)
        
        
        
        #adding enumeration
        water_bonds2 = adding_numeration_to_items(water_bonds)
        water_angles2 = adding_numeration_to_items(water_angles)
        num_bond = water_bonds2.show_last_number()
        num_angle = water_angles2.show_last_number()
        
    
        water_bonds_ready = water_bonds2.add_enumeration()
        water_angles_ready = water_angles2.add_enumeration()

            
        water_bonds_types = water_bonds2.count_bond_types()
        water_angles_types = water_angles2.count_bond_types()
        
        #DATA FILE:
        x = data_file(num_atoms, num_bond,  num_angle, atom_types_dic, water_bonds_types, water_angles_types, dimen, atom_list,
                      water_bonds_ready, water_angles_ready)
        y = x.produce_data_file()
        for f in y:
            for m in f:
                print(*m,  file = outfile)
                
        
        #anchors in nanotube
        atoms_in_anchors = getting_anchors_evolution(atom_list, nanchors, dimen, atoms_per_anchor)
        print(*atoms_in_anchors, file = anchors_sicnt_T)
#         for ids_in in atoms_in_anchors:
#             print(ids_in, file = anchors_cnt_T)
            #print(*atoms_in_anchors, file = anchors_T)
            
        sicnt_wo_anchors = sicnt_no_anchors(sicnt_s, atoms_in_anchors)
        print(*sicnt_wo_anchors, file = sicar_no_anch_out)
        #print(*cnt_no_anchors, file = car_no_anch_out)
        
        almost_all = sicnt_wo_anchors + nofrozm #water_s
        print(*almost_all, file = sicnt_no_anchors_and_water_out)
        
        
        splitted_anchors_group = split_list(atoms_in_anchors)
        #print(*splitted_anchors_group[0])
        #print(*splitted_anchors_group[1])
        #2 groups with half of the anchors;
        #there are 3 anchors, so each group will 'have' 1.5 anchor
        
        sicnt_no_anchor1 = sicnt_no_anchors(sicnt_s, splitted_anchors_group[0])
        sicnt_no_anchor2 = sicnt_no_anchors(sicnt_s, splitted_anchors_group[1])
        
        print(*sicnt_no_anchor1, file = sicnt_group1) #atoms from anchors 1 removed
        print(*sicnt_no_anchor2, file = sicnt_group2) #atoms from anchors 2 removed
        
        
        
#####################################################################
#--------------------FILES FOR INPUT FILE---------------------------#


def split_list(a_list):
    
    splited_in2 = np.array_split(a_list, 2)
    return splited_in2
    
    
def file_to_list(infile):
    out_list = []
    for lines in infile:
        lines_s = int(lines.strip())
        out_list.append(lines_s)
    return out_list

def sicnt_no_anchors(master_list, anchor_list):
    
    master_list1 = set(master_list)#.read().split()
    anchor_list1 = set(anchor_list)#.read().split()
    
    file1_minus_file2 = list(set(master_list1) - set(anchor_list1))

    #file1_minus_file2 = list(set(master_list) - set(anchor_list))
    
    return file1_minus_file2

def get_ids_from_dist_lists(array_l):
    final_anchors = []
    for i in array_l:
        final_anchors.append(i[1])
        
    return final_anchors
        

def getting_anchors_evolution(atom_list, number_of_anchors, dimensions, atoms_per_anchor):
    #number_of_anchors_desired, atom_list, atom_type1, atom_type2
   
    test_dic = {}
    xyx_sicnt = []
    anchors = []
    si_anchor = []
    c_anchor = []
    z_dim = dimensions[2]
    bond_length = 1.80
    
    #return a list
    for lines_of_atoms in atom_list:
        if lines_of_atoms[2] == 1 or lines_of_atoms[2] == 2:
            to_append = [lines_of_atoms[0], lines_of_atoms[2], lines_of_atoms[4], lines_of_atoms[5], lines_of_atoms[6]]
            xyx_sicnt.append(to_append)
    #return xyx_sicnt
    min_z = min(inner_list[4] for inner_list in xyx_sicnt)
    max_z = max(inner_list[4] for inner_list in xyx_sicnt)
    
    range_of_z = (max_z - min_z)/float(number_of_anchors+1)
    
    multi_list = [] 
    #short list with the 3 positions in z for the anchors. The atoms chosen must be the nearest to these points
    for x in range(1, number_of_anchors+1, 1):
        multiplier = x*range_of_z
        multi_list.append(multiplier)
    #print(atom_list)
    
#     if number_of_anchors == 3:
#         z_closest1 = min(xyx_sicnt, key=lambda z:abs(z[4]-multi_list[0]))
#         z_closest2 = min(xyx_sicnt, key=lambda z:abs(z[4]-multi_list[1]))
#         z_closest3 = min(xyx_sicnt, key=lambda z:abs(z[4]-multi_list[2]))
#         #return z_closest
# 
#     for lines_of_atoms in atom_list:
#         if lines_of_atoms[6] == z_closest1[4] or lines_of_atoms[6] == z_closest2[4] or \
#             lines_of_atoms[6] == z_closest3[4] or lines_of_atoms[6] == (z_closest1[4] - bond_length) or \
#             lines_of_atoms[6] == (z_closest2[4] - bond_length) or lines_of_atoms[6] == (z_closest3[4] - bond_length):
#             
#             if lines_of_atoms[2] == 1:
#                 si_anchor.append(lines_of_atoms[0])
#             elif lines_of_atoms[2] == 2:
#                 c_anchor.append(lines_of_atoms[0])
#     anchor_ids = c_anchor + si_anchor
#     anchor_ids_s = sorted(anchor_ids)
    
    
    distances_for_anchors_1_C, distances_for_anchors_1_Si, distances_for_anchors_2_C, distances_for_anchors_2_Si, \
    distances_for_anchors_3_C, distances_for_anchors_3_Si =  ([] for i in range(6))
    
  
    if number_of_anchors == 3:
        for i in atom_list: 
            dist_1 = abs(i[6] - multi_list[0])
            dist_2 = abs(i[6] - multi_list[1])
            dist_3 = abs(i[6] - multi_list[2])
            if i[2] == 1: #type 1 Si
                to_append_1_Si = [dist_1 , i[0]] #distance, id
                distances_for_anchors_1_Si.append(to_append_1_Si)
                
                to_append_2_Si = [dist_2 , i[0]] #distance, id
                distances_for_anchors_2_Si.append(to_append_2_Si)
                
                to_append_3_Si = [dist_3 , i[0]] #distance, id
                distances_for_anchors_3_Si.append(to_append_3_Si)
                
                
            elif i[2] == 2: #type 2 C
                to_append_1_C = [dist_1 , i[0]]
                distances_for_anchors_1_C.append(to_append_1_C)
                
                to_append_2_C = [dist_2 , i[0]]
                distances_for_anchors_2_C.append(to_append_2_C)
                
                to_append_3_C = [dist_3 , i[0]]
                distances_for_anchors_3_C.append(to_append_3_C)
                
                
        distances_for_anchors_1_Si_s = sorted(distances_for_anchors_1_Si)
        distances_for_anchors_1_C_s = sorted(distances_for_anchors_1_C)
        
        distances_for_anchors_2_Si_s = sorted(distances_for_anchors_2_Si)
        distances_for_anchors_2_C_s = sorted(distances_for_anchors_2_C)
        
        distances_for_anchors_3_Si_s = sorted(distances_for_anchors_3_Si)
        distances_for_anchors_3_C_s = sorted(distances_for_anchors_3_C)
#------------------
        
        si_anchor_1_ids = get_ids_from_dist_lists(distances_for_anchors_1_Si_s)
        c_anchor_1_ids =  get_ids_from_dist_lists(distances_for_anchors_1_C_s)
        
        si_anchor_2_ids = get_ids_from_dist_lists(distances_for_anchors_2_Si_s)
        c_anchor_2_ids =  get_ids_from_dist_lists(distances_for_anchors_2_C_s)
        
        si_anchor_3_ids = get_ids_from_dist_lists(distances_for_anchors_3_Si_s)
        c_anchor_3_ids =  get_ids_from_dist_lists(distances_for_anchors_3_C_s)
        
        
        
    
    #print(distances_for_anchors_1_Si_s)
    
    anchor1 = si_anchor_1_ids[:int(atoms_per_anchor/2)] + c_anchor_1_ids[:int(atoms_per_anchor/2)]
    anchor2 = si_anchor_2_ids[:int(atoms_per_anchor/2)] + c_anchor_2_ids[:int(atoms_per_anchor/2)]
    anchor3 = si_anchor_3_ids[:int(atoms_per_anchor/2)] + c_anchor_3_ids[:int(atoms_per_anchor/2)]
            
    anchor_ids = anchor1 + anchor2 + anchor3
    anchor_ids_s = sorted(anchor_ids)
            
    
            
            
         
    return anchor_ids_s


#####################################################################
#-----------------------------DATA-FILE-----------------------------#
class data_file(object):
    
    def __init__(self, num_atoms, num_bonds, num_angles, 
                 num_atom_types, num_bond_types, num_angle_types, 
                 dimensions, atom_list, water_bonds_ready, water_angles_ready, num_dihedrals = 0, num_impropers = 0):
        
        self._num_atoms = num_atoms
        self._num_bonds = num_bonds
        self._num_angles = num_angles
        self._num_dihedrals = num_dihedrals
        self._num_impropers = num_impropers
        
        self._num_atom_types = num_atom_types
        self._num_bond_types = num_bond_types
        self._num_angle_types = num_angle_types
        
        self._water_bonds_ready = water_bonds_ready
        self._water_angles_ready = water_angles_ready
        
        self._atom_list = atom_list
        self._dimensions = []
        
        self._atoms = 'atoms'
        self._angles = 'angles' 
        self._bonds = 'bonds'
        self._dihedrals = 'dihedrals'
        self._impropers = 'impropers'
        
        self._atom_types = 'atom types'
        self._bond_types = 'bond types'
        self._angle_types = 'angle types'
        #not used here:
        self._dihedral_types = 'dihedral types'
        self._improper_types = 'improper types'
        
        self._x = dimensions[0]
        self._y = dimensions[1]
        self._z = dimensions[2]
        
        self._x0 = 0.0
        self._y0 = 0.0
        self._z0 = 0.0
        
        self._x_tag =   'xlo xhi'
        self._y_tag =   'ylo yhi'
        self._z_tag =   'zlo zhi'
        
        #self._masses = {'1': 28.0850 '2': 12.011, '3': 1.008, '4':15.9994}
        #self._masses = 'Masses'
        self._mass1 = 1, 28.0850
        self._mass2 = 2, 12.011
        self._mass3 = 3, 1.008
        self._mass4 = 4, 15.994, "\n"
        
    def produce_data_file(self):
                    
        h0 = [
        ['This is a datafile \n'],
        [self._num_atoms , self._atoms],
        [self._num_bonds , self._bonds],
        [self._num_angles , self._angles],
        [self._num_dihedrals , self._dihedrals],
        [self._num_impropers , self._impropers, "\n"],
        
        [self._num_atom_types, self._atom_types],
        [self._num_bond_types, self._bond_types],
        [self._num_angle_types, self._angle_types, "\n"],
        
        [self._x0, self._x, self._x_tag],
        [self._y0, self._y, self._y_tag],
        [self._z0, self._z, self._z_tag, "\n"]]
        
        
        body_data_file = [h0, [['Masses \n']], [self._mass1, self._mass2, self._mass3, self._mass4], [['Atoms \n']], self._atom_list, ' ', [['Bonds \n']], 
        self._water_bonds_ready, ' ', [['Angles \n']], self._water_angles_ready] 
        return body_data_file
        #return self._water_bonds_ready
    
#####################################################################
#-----------CLASESS AND METHODS FOR BONDS AND ANGLES----------------#
            
class getting_water_items(object):
    
    def __init__(self, bond_list, angle_list):
        self._bond_list = bond_list
        self._angle_list = angle_list
        
        self._water_bonds = []
        self._water_angles = []
        
        self._bond_water_type = int(1)
        self._angle_water_type = int(1)
        
    def get_water_bonds(self):
        for bonds_pairs in self._bond_list:
            if bonds_pairs[0] == self._bond_water_type:
                 self._water_bonds.append(bonds_pairs)
        return self._water_bonds
    
    def get_water_angles(self):
        for angles in self._angle_list:
            if angles[0] == self._angle_water_type:
                 self._water_angles.append(angles)
        return self._water_angles
    
    
#----------------------------------------------------------------#
class adding_numeration_to_items(object):
    
    def __init__(self, list_of_items):
        self._returned_list = []
        self._input_list = list_of_items
        
    def add_enumeration(self):
    
        counter_of_items = 0
        for lists in self._input_list:
            if lists[0] == 1:
                counter_of_items += 1
                lists.insert(0, counter_of_items)
                self._returned_list.append(lists)
                #print(lists)
        return self._returned_list
    
    def show_last_number(self):
        
        counter_of_items = 0
        for lists in self._input_list:
                counter_of_items += 1
        return counter_of_items
    
    def count_atom_types(self):
        counter_of_items = 0
        a_list = []
        for lists in self._input_list:
                a_list.append(lists[2])
        res = Counter(a_list)
        #return res
        for key, value in res.items():
            counter_of_items += 1
        return counter_of_items
    
    def count_bond_types(self):
        counter_of_items = 0
        a_list = []
        for lists in self._input_list:
                a_list.append(lists[1])
        res = Counter(a_list)
        #return res
        for key, value in res.items():
            counter_of_items += 1
        return counter_of_items
        
    
    def count_angle_types(self):
        counter_of_items = 0
        a_list = []
        for lists in self._input_list:
                a_list.append(lists[1])
        res = Counter(a_list)
        #return res
        for key, value in res.items():
            counter_of_items += 1
        return counter_of_items
        
#----------------------------------------------------------#

##########################################################
#---------------------BONDS------------------------------#
            
def getting_bonds(conections, carbon, silicon, hydrogen, oxygen):
            final_bond_list = []
            bonds_list = []
            #getting bonds:
            for conecti in conections:
                bond_object = bond()  #this is the instance of an object, a method will create the object
                
                bond_created = bond_object.create_bond(conecti, carbon, hydrogen, oxygen)
                if bond_created != None:
                    bx = tuple(bond_created)
                    bonds_list.append(bx)

            #eliminating duplicate bonds:
            bonds_no_duplicates = sorted(set(bonds_list))
    
            #print(bonds_no_duplicates)
            bond_type1 = int(1) #O-H
            bond_type2 = int(2) #Si-C
            for bond_tup in bonds_no_duplicates:
                bonds_pairs = list(bond_tup)
                #print(bonds_list_no_dup)
                
                if bonds_pairs[0] in carbon:
                    bonds_pairs.insert(0,bond_type2)
                elif bonds_pairs[0] in oxygen:
                    bonds_pairs.insert(0,bond_type1)
                elif bonds_pairs[0] in silicon:
                    bonds_pairs.insert(0,bond_type2)
                final_bond_list.append(bonds_pairs)
            #print(final_bond_list)
            return final_bond_list
            
    
class bond(object):
    #an angle is formed by 3 atoms
    #angles types: water: 1, cnt: 2
    def __init__(self):
        
        self._single_bonds = []
        self._bond_type1 = int(1)
        self._bond_type2 = int(2)
    
    def create_bond(self, items_raw_bonds, carbon, hydrogen, oxygen):
        #cou = 0
        if len(items_raw_bonds) == 3:
            items_for_bonds_0 = int(items_raw_bonds[0])
            items_for_bonds_1 = int(items_raw_bonds[1])
            items_for_bonds_2 = int(items_raw_bonds[2])
            
            #add the type of angle
            if items_for_bonds_0 in oxygen:
                single_bond1 = sorted([items_for_bonds_0, items_for_bonds_1])
                single_bond2 = sorted([items_for_bonds_0, items_for_bonds_2])
                self._single_bonds =  [single_bond1, single_bond2]
                #return self._single_bonds
            
            elif items_for_bonds_0 in carbon:
                single_bond1 = sorted([items_for_bonds_0, items_for_bonds_1])
                single_bond2 = sorted([items_for_bonds_0, items_for_bonds_2])
                self._single_bonds =  [single_bond1, single_bond2]
                #return self._single_bonds
            
            for items_in in self._single_bonds:
                return items_in

            
        elif len(items_raw_bonds) == 4:
            items_for_bonds_0 = int(items_raw_bonds[0])
            items_for_bonds_1 = int(items_raw_bonds[1])
            items_for_bonds_2 = int(items_raw_bonds[2])
            items_for_bonds_3 = int(items_raw_bonds[3])
            
            if items_for_bonds_0 in oxygen:
                single_bond1 = sorted([items_for_bonds_0, items_for_bonds_1])
                single_bond2 = sorted([items_for_bonds_0, items_for_bonds_2])
                single_bond3 = sorted([items_for_bonds_0, items_for_bonds_3])
                self._single_bonds =  [single_bond1, single_bond2, single_bond3]
                #return self._single_bonds
            
            elif items_for_bonds_0 in carbon:
                single_bond1 = sorted([items_for_bonds_0, items_for_bonds_1])
                single_bond2 = sorted([items_for_bonds_0, items_for_bonds_2])
                single_bond3 = sorted([items_for_bonds_0, items_for_bonds_3])
                self._single_bonds =  [single_bond1, single_bond2]
                #return self._single_bonds
            
            for items_in in self._single_bonds:
                return items_in
            
        elif len(items_raw_bonds) == 2:
            items_for_bonds_0 = int(items_raw_bonds[0])
            items_for_bonds_1 = int(items_raw_bonds[1])
            
            if items_for_bonds_0 in oxygen or items_for_bonds_0 in hydrogen:
                single_bond = sorted([items_for_bonds_0, items_for_bonds_1])
                return single_bond
            
            elif items_for_bonds_0 in carbon:
                single_bond = sorted([items_for_bonds_0, items_for_bonds_1])
                return single_bond
            

###########################################################
#--------------------ANGLES-------------------------------#

def getting_angles(conections, carbon, hydrogen, oxygen):
            angles_list = []
            #getting the angles:
            for conecti in conections:
                angle_object = angle()
                angle_created = angle_object.create_angle(conecti, carbon, hydrogen, oxygen)
                if angle_created != None:
                    angles_list.append(angle_created)
                    #print(angle_created)
            return angles_list
        
#----------------------------------------------------------#

class angle(object):
    #an angle is formed by 3 atoms
    #angles types: water: 1, cnt: 2
    def __init__(self):
        
        self._single_angles = []
        self._ang_type1 = int(1)
        self._ang_type2 = int(2)
    
    def create_angle(self, items_raw_angles, carbon, hydrogen, oxygen):
        #cou = 0
        if len(items_raw_angles) == 3:
            items_for_angles_0 = int(items_raw_angles[0])
            items_for_angles_1 = int(items_raw_angles[1])
            items_for_angles_2 = int(items_raw_angles[2])
            
            #add the type of angle
            if items_for_angles_0 in oxygen:
                single_angle = [self._ang_type1, items_for_angles_1, items_for_angles_0, items_for_angles_2] 
                return single_angle
                
            
            elif items_for_angles_0 in carbon:
                single_angle = [self._ang_type2, items_for_angles_1, items_for_angles_0, items_for_angles_2] 
                return single_angle
            
        elif len(items_raw_angles) == 4:
            items_for_angles_0 = int(items_raw_angles[0])
            items_for_angles_1 = int(items_raw_angles[1])
            items_for_angles_2 = int(items_raw_angles[2])
            items_for_angles_3 = int(items_raw_angles[3])
            
            if items_for_angles_0 in oxygen:
                single_angle1 = [self._ang_type1, items_for_angles_1, items_for_angles_0, items_for_angles_2] 
                single_angle2 = [self._ang_type1, items_for_angles_1, items_for_angles_0, items_for_angles_3]
                single_angle3 = [self._ang_type1, items_for_angles_3, items_for_angles_0, items_for_angles_2]
                self._single_angles = [single_angle1, single_angle2, single_angle3] 
            
            elif items_for_angles_0 in carbon:
                single_angle1 = [self._ang_type2, items_for_angles_1, items_for_angles_0, items_for_angles_2] 
                single_angle2 = [self._ang_type2, items_for_angles_1, items_for_angles_0, items_for_angles_3]
                single_angle3 = [self._ang_type2, items_for_angles_3, items_for_angles_0, items_for_angles_2]
                self._single_angles = [single_angle1, single_angle2, single_angle3] 
            
            for items_in in self._single_angles:
                return items_in

###########################################################
#------------------ATOMIC-POSITIONS----------------------#
class split_types(object):
    def __init__(self,id, element):
        self._elem = element
        self._id = int(id)
        
    def get_Si(self):
        if self._elem == 'Si':
            return self._id
        
    def get_C(self):
        if self._elem == 'C':
            return self._id
        
    def get_SiC(self):
        if self._elem == 'C' or self._elem == 'Si':
            return self._id
            
    def get_H(self):  
        if self._elem == 'H':
            return self._id
            #print(hydrogen)
    def get_O(self):
        if self._elem == 'O':
            return self._id
    
    def count_atoms(self):
        counter_of_items = 0
        for lists in self._input_list:
            if lists[0] == 1:
                counter_of_items += 1
        return counter_of_items
    
    def masses(self):
        
        
        pass
#---------------------------------------------------------#            
             
class atom_line(object):
    
    def __init__(self, id, x, y, z, element):
        self._id = int(id)
        self._x = float(x) 
        self._y = float(y)
        self._z = float(z) 
        self._elem = element 
        self._line_list = []
        self._fin_list = []
        
    #CHARGES FOR TIP4P-EW MODEL
    def get_q_type_mol_tag(self):
        #atoms types: Si 1, C 2, H 3, O 4:
        #molecule tag: 1 SiCNT, 2 water
        #some important parameters like LJ must be written in the Input file
        #Parameters based in TIP4P/ice
        if self._elem == 'Si':
            self._q = 0.60
            self._type = int(1)
            self._mole = int(1)
            line_list = [self._id, self._mole, self._type, self._q, self._x, self._y, self._z]
            return line_list
         
        elif self._elem == 'C':
            self._q = -0.60
            self._type = int(2)
            self._mole = int(1)
            line_list = [self._id, self._mole, self._type, self._q, self._x, self._y, self._z]
            return line_list
             
        elif self._elem == 'H':
            self._q = 0.5897
            self._type = int(3)
            self._mole = int(2)
            line_list = [self._id, self._mole, self._type, self._q, self._x, self._y, self._z]
            return line_list
        
        elif self._elem == 'O':
            self._q = -1.1794
            self._type = int(4)
            self._mole = int(2)
            line_list = [self._id, self._mole, self._type, self._q, self._x, self._y, self._z]
            return line_list

if __name__ == "__main__": main()
