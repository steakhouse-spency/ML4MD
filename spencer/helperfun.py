import numpy as np
import math as mt
import random
import glob
import os
import copy

from globalvar import *

def seperate_empty_tube(empty_tube_file):
    header = None
    atoms = []
    conect = []

    # get first line as header
    header = empty_tube_file.readline()

    # iterate thru file, seperate according to first column
    for i, row in enumerate(empty_tube_file):
        col = row.split()

        if col[0] == residue:
            # handle when col 4 and col 5 get combined
            print(col[4])

            if len(col[4]) > 1:
                atoms.append(col[4][0])
                atoms.append(col[4][1:])
                print( col[4][1:])
            else:
                atoms.append(col)
        elif col[0] == "CONECT":
            conect.append(col)
        elif col[0] == "END":
            break
        else:
            print("file format error: check empty tube file")
            exit(1)

    exit(1)


    # move atoms/conect coordinates
    xyz=(100,100,420)

    # find min (x, y, z) only if negative
        # add min to all atoms

    normalizeCoord(atoms)

    # 
    # min(min(atoms))


    return header, atoms, conect


def normalizeCoord(atoms):
    x_min = 10000
    y_min = 10000
    z_min = 10000
    for atom in atoms:
        x_min = min(x_min, float(atom[6]))
        y_min = min(y_min, float(atom[7]))
        z_min = min(z_min, float(atom[8]))

    if x_min > 0: x_min = 0
    if y_min > 0: y_min = 0
    if z_min > 0: z_min = 0

    x_min *= -1
    y_min *= -1
    z_min *= -1

    for atom in atoms:
        atom[6] = float(atom[6]) + x_min
        atom[7] = float(atom[7]) + y_min
        atom[8] = float(atom[8]) + z_min





def assign_positions_to_oxygens_list_pdb_new(number_of_atoms, tube_center, startn):
    
    chuncks_z = L/side_cube
    N_z = number_of_atoms/chuncks_z
    #print('N_z = ', N_z)
    N_z_f = mt.modf(N_z)
    #0 is decimal and 1 is the integer
    # print(N_z_f[0], N_z_f[1], chuncks_z)
    
    every_ring = int(1/N_z_f[0])
    atom_per_ring = int(N_z_f[1]) #ideal situation
    

    # ASK JOSE: do we want L updated??

    list_positions_0 = []
    list_positions_final = []
    
    r= (D/2.0)- safety_clearance #safety clearance
    
    
    count_additional = 1
    zin = 1.0
    
    #pentagon
    while zin <= L:
        for n in range(0,5):
            x = r*0.6*mt.cos(mt.radians(90+n*72))
            y = r*0.6*mt.sin(mt.radians(90+n*72))

            id = startn
            #increase_by = 1.0
            new_x = x + tube_center[0]
            new_y = y + tube_center[1]
            new_z = zin
            
            coords_0 = [id, [new_x,new_y, new_z]]

            startn = startn + 3
            
            list_positions_0.append(coords_0)
            #zin = zin + side_cube
            #
        zin = zin + side_cube
        count_additional = count_additional + 1

        if count_additional%4 == 0: #because of the 0.25  decimal
            zin = zin - side_cube
            id = startn

            new_x = tube_center[0]
            new_y = tube_center[1]
            new_z = zin
            coords_0 = [id, [new_x,new_y, new_z]]
            
            list_positions_0.append(coords_0)
            #zin = zin + side_cube
            #count_additional = count_additional + 1
            startn = startn + 3
            zin = zin + side_cube

    for i in list_positions_0:
        list_positions_final.append(i[1])
    
    return list_positions_final
        
def lj_oxy(r, sigma_O, eps_O):
        pe_lj = 4*(eps_O)* ((pow((sigma_O/r),12)) - (pow((sigma_O/r),6)))
        return pe_lj
        
def pe_pairwise_LJ(list_oxy_positions, sigma_O, eps_O, cutoff):
    #interactions between O and H = 0
    list_b = copy.deepcopy(list_oxy_positions)
    eps_mixed = 1
    sig_mixed = 1
    #pe_lj = 4*(sig_mixed)* ((pow(sig_mixed/r),12) - (pow(sig_mixed/r),6))
    
    n = len(list_oxy_positions)
    pe_sum = 0.0
    for i in list_oxy_positions:
        for j in list_b:
            if i[0] != j[0]:
                r_test = (j[1][0] - i[1][0])*(j[1][0] - i[1][0]) + (j[1][1] - i[1][1])*(j[1][1] - i[1][1]) + \
                             (j[1][2] - i[1][2])*(j[1][2] - i[1][2])
                             
                if r_test <= cutoff*cutoff:
                    r = np.sqrt(r_test)
                    pe_ij = lj_oxy(r, sigma_O, eps_O)
                    
                    pe_sum += pe_ij
                    
                    
    return pe_sum
                    

# only needs 
def get_tube_center(tube_atoms):

        coordinatesxy = []

        for i, atom in enumerate(tube_atoms):
            # 6 = x, 7 = y, 8 = z
            coordsxy = [float(atom[6]), float(atom[7]), float(atom[8])]
            coordinatesxy.append(coordsxy)


# ASK JOSE:  coordinates need to be scaled by this point
       
        # randomly choose 3 (x,y,z) coordinates 
        xyz1 = random.choice(coordinatesxy)
        xyz2 = random.choice(coordinatesxy)
        xyz3 = random.choice(coordinatesxy)

        # only compare x and y?
        while xyz1[0] == xyz2[0] or xyz1[0] == xyz3[0] or xyz2[0] == xyz3[0] or \
              xyz1[1] == xyz2[1] or xyz1[1] == xyz3[1] or xyz2[1] == xyz3[1]:
            
            xyz1 = random.choice(coordinatesxy)
            xyz2 = random.choice(coordinatesxy)
            xyz3 = random.choice(coordinatesxy)
        
        circle_points_3 = [xyz1, xyz2, xyz3]
        circle_points_3_s = sorted(circle_points_3, key = lambda list_of_3:list_of_3[0])
        #print(circle_points_3_s)
        
        # geometry to compute center
        p1 = circle_points_3_s[0]
        p2 = circle_points_3_s[1]
        p3 = circle_points_3_s[2]
        ma = (p2[1] - p1[1]) / (p2[0] - p1[0])
        mb = (p3[1] - p2[1]) / (p3[0] - p2[0])
        center_x = ((ma*mb*(p1[1] - p3[1])) + (mb*(p1[0] + p2[0])) - (ma*(p2[0] + p3[0]))) /(2*(mb-ma))
        center_y = ((-1/ma) * (center_x -((p1[0] + p2[0])/2.0))) + ((p1[1] + p2[1])/2.0)
        
        return [center_x, center_y]
        
def water_conections(oxygen_quantity, startn):
    conect_list = []
    i = startn #1st atom, oxygen
    conect = 'CONECT'
    while i < ((3*oxygen_quantity)+startn):
        conection1 = [conect, i, i+1, i+2]
        conection2 = [conect, i+1, i]
        conection3 = [conect, i+2, i]
        cnt = [conection1, conection2, conection3]
        i = i + 3 #to skip the previous atoms
        #print(i)
        conect_list.extend(cnt)
    #print(oxygen_quantity+startn)
    return conect_list
    
    
class water_class_pdb(object):
    def __init__(self, oxygen_list, typeO, hydrogen1, hydrogen2, typeH, start_number):
        self._Otype = typeO
        self._Htype = typeH
        self._oxy_list = oxygen_list
        self._h1_list = hydrogen1
        self._h2_list = hydrogen2
        self._start = int(start_number)
        self._atoms = residue
        self._xxx = 'XXX'
        self._anumb = 'MOL'
        self._one = '1.00'
        self._zero = '0.00'
        self._oxy_element = 'O'
        self._hyd_element = 'H'
        
    def oxy_processed(self):
        oxy_final = []
        n = self._start
        for oxy_coords in self._oxy_list:
            final_line = [self._atoms, n, self._anumb, self._xxx, n, format(oxy_coords[0], '.3f'),\
                           round(oxy_coords[1],3), round(oxy_coords[2],3),  self._one,  self._zero,  self._oxy_element]
            n = n+3
            oxy_final.append(final_line)
        return oxy_final
    
    def hydro_processed1(self):
        hydro_final1 = []
        nh1 = self._start + 1
        for hydro1_coords in self._h1_list:
            final_line = [self._atoms, nh1, self._anumb, self._xxx, nh1-1, format(hydro1_coords[0],'.3f'),\
                           format(hydro1_coords[1],'.3f'), format(hydro1_coords[2],'.3f'),  self._one,  self._zero,  self._hyd_element]
            nh1 = nh1 + 3
            hydro_final1.append(final_line)
        return hydro_final1
    
    def hydro_processed2(self):
        hydro_final2 = []
        nh2 = self._start + 2
        for hydro2_coords in self._h2_list:
            final_line = [self._atoms, nh2, self._anumb, self._xxx, nh2-2, format(hydro2_coords[0],'.3f'),\
                           format(hydro2_coords[1],'.3f'), format(hydro2_coords[2],'.3f'),  self._one,  self._zero,  self._hyd_element]
            nh2 = nh2 + 3
            hydro_final2.append(final_line)
        return hydro_final2

    
def final_position_for_H(oxy_list, a_list):
    final_position_of_a = []
    for a1, a0 in zip(oxy_list, a_list):
        x = a1[0] + a0[0]
        y = a1[1] + a0[1]
        z = a1[2] + a0[2]
        coords = [x, y, z]
        final_position_of_a.append(coords)
        
    return final_position_of_a

def calculate_water_molecules_inside_nanotube(density_g_cc):
    water_mass = 18.01528
    avogadro = 6.022e23
    cm3_ti_ang3 = 1e24
    R = D/2.0
    v_ang3 = R*R*L*np.pi
    den_molecules_ang3 =    avogadro*density_g_cc/(water_mass*cm3_ti_ang3)
    # print('vang3 = ', v_ang3)
    # print(den_molecules_ang3)
    water_molecules_to_fill = den_molecules_ang3*v_ang3
    return int(water_molecules_to_fill)



def adding_H1(number_of_atoms, angle, bond_length):
    rand_rot = np.random.randint(1, 3)
    values = [0,0,0]
    totalh1 = []
    #sprint(values)
    for _numbers_atoms in range((1*number_of_atoms)+1):
        rand_ini = np.random.randint(1, 7)
        #print(rand_ini)
        if rand_ini == 1 and rand_rot == 1:
            basic_vector = [0.0, 0.0, bond_length]
        
        elif rand_ini == 1 and rand_rot == 2:
            basic_vector = [0.0, 0.0, bond_length]
        #===========================================================================
        elif rand_ini == 2 and rand_rot == 1:
            basic_vector = [0.0, bond_length, 0.0]
            
        elif rand_ini == 2 and rand_rot == 2:
            basic_vector = [0.0, bond_length, 0.0]

        #===============================================================================
        elif rand_ini == 3 and rand_rot == 1:
            basic_vector = [bond_length, 0.0, 0.0]
            
        elif rand_ini == 3 and rand_rot == 2:
            basic_vector = [bond_length, 0.0, 0.0]
                
        ##########################=====================================================
        elif rand_ini == 4 and rand_rot == 1:
            basic_vector = [0.0, 0.0, -bond_length]
            
        elif rand_ini == 4 and rand_rot == 2:
            basic_vector = [0.0, 0.0, -bond_length]
            
        #===========================================================================
        elif rand_ini == 5 and rand_rot == 1:
            basic_vector = [0.0, -bond_length, 0.0]
            
        elif rand_ini == 5 and rand_rot == 2:
            basic_vector = [0.0, -bond_length, 0.0]
            
        #===============================================================================
        elif rand_ini == 6 and rand_rot == 1:
            basic_vector = [-bond_length, 0.0, 0.0]
            
        elif rand_ini == 6 and rand_rot == 2:
            basic_vector = [-bond_length, 0.0, 0.0]
             
        totalh1.append(basic_vector)
    return totalh1

def adding_H2(h1_list, angle, bond_length):
    r_angle = mt.radians(angle)
    cos_ang = np.cos(r_angle)
    total_h2 = []

    for basic_vector in h1_list:
        #print(basic_vector)
        
        rand_rot = np.random.randint(1, 3)

        if basic_vector[0] == 0.0 and basic_vector[1] == 0.0:
            #print(basic_vector)
            if rand_rot == 1:
                xh2 = basic_vector[0]
                yh2 = (basic_vector[1]*np.cos(r_angle) - basic_vector[2]*np.sin(r_angle))
                zh2 = (basic_vector[1]*np.sin(r_angle) + basic_vector[2]*np.cos(r_angle))
            elif rand_rot == 2:
                xh2 = basic_vector[0]*np.cos(r_angle) + basic_vector[2]*np.sin(r_angle)
                yh2 = basic_vector[1]
                zh2 = -basic_vector[0]*np.sin(r_angle) + basic_vector[2]*np.cos(r_angle) 
                
        
        elif basic_vector[0] == 0.0 and basic_vector[2] == 0.0:
            if rand_rot == 1:
                xh2 = basic_vector[0]
                yh2 = (basic_vector[1]*np.cos(r_angle) - basic_vector[2]*np.sin(r_angle))
                zh2 = (basic_vector[1]*np.sin(r_angle) + basic_vector[2]*np.cos(r_angle))
            elif rand_rot == 2:
                xh2 = basic_vector[0]*np.cos(r_angle) - basic_vector[1]*np.sin(r_angle)
                yh2 = basic_vector[0]*np.sin(r_angle) + basic_vector[1]*np.cos(r_angle)
                zh2 = basic_vector[2]
                
        
        elif basic_vector[1] == 0.0 and basic_vector[2] == 0.0:
            if rand_rot == 1:
                xh2 = basic_vector[0]*np.cos(r_angle) + basic_vector[2]*np.sin(r_angle)
                yh2 = basic_vector[1]
                zh2 = -basic_vector[0]*np.sin(r_angle) + basic_vector[2]*np.cos(r_angle) 
                
            elif rand_rot == 2:
                xh2 = basic_vector[0]*np.cos(r_angle) - basic_vector[1]*np.sin(r_angle)
                yh2 = basic_vector[0]*np.sin(r_angle) + basic_vector[1]*np.cos(r_angle)
                zh2 = basic_vector[2]
                
                
        vector_h2 = [xh2, yh2, zh2]
        total_h2.append(vector_h2)
        
        
    return total_h2