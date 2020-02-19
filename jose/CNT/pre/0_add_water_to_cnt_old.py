
#import fileinput
import numpy as np
import math as mt
import random
import glob
import os
def main():
    # This file fill an empty nanotube with water. The input file required is only a pdb file with the coordinates of the nanotube
    # and the 'CONECT' section.
    
    #density:
    rho = 1.0
    
    #nanotube dimensions
    D = 10.6
    L = 288.378
    
    #water parameters:
    bond_length = 0.9572
    angle = 104.52
    safety_clearance = 1.5
    
    
    with open('water_pdb.txt', 'a+') as waterout, open('water_conect.txt', 'w+') as water_conect_out, \
         open('tube_pdb.txt', 'w+') as tubeout, open('tube_conect.txt', 'w+') as tube_conect_out, \
         open('header_section.txt', 'w+') as header_out, \
        open('150xy-long-10.6.pdb', 'r') as file_get_center_tube, open('final_pdb_{}.pdb'.format(rho), 'w') as final_output:
        #the file file_get_center_tube is just a nanotube, nothing else, in pdb format
        
        #file_get_center_tube = open('pdb_sicnt_prueba.pdb', 'r')
        xy_center_tube = gettin_center_of_tube(file_get_center_tube)
        file_get_center_tube.seek(0)
#getting the pieces from the pdb file:
        file_obj = get_list_for_pdb_input_file(file_get_center_tube)
        header_part = file_obj.getting_header_from_pdb_file()
        atoms_tube = file_obj.getting_Atoms()
        conect_tube = file_obj.getting_conect()
        
        for lines in header_part:
            print(*lines, file = header_out)
            
        for lines in conect_tube:
            print(*lines, file = tube_conect_out)
            
        for j in atoms_tube:
            new_line = '{:<0} {:>6} {:>4} {:>3} {:>5} {:>11} {:>7} {:>7} {:>5} {:>5} {:>11}'.format( \
                        j[0], j[1], j[2], j[3], j[4], j[5], j[6], j[7], j[8], j[9], j[10])
            print(new_line, file = tubeout)
            
#getting the water. 
        cylinder_x = xy_center_tube[0] #center of the cylinder
        cylinder_y = xy_center_tube[1]
        startn = 1 + len(atoms_tube)
        #2881
        oxygen_quantity = int(calculate_water_molecules_inside_nanotube(D, L, rho))
        print(oxygen_quantity)
        oxygens = assign_positions_to_oxygens_list_pdb(D, L, oxygen_quantity, bond_length, cylinder_x, cylinder_y, safety_clearance)
            
        hydro1 = adding_H1(oxygen_quantity, angle, bond_length)
        #print(hydro1)
        hydro2 = adding_H2(hydro1, angle, bond_length)
        #print(hydro2)
        #print(*hydro, file = fileout2)
        #print(hydro)
        h1_final = final_position_for_H(oxygens, hydro1) #move the water atoms from the origin to the random positions generated in oxygens
        h2_final = final_position_for_H(oxygens, hydro2)

        water_positions = water_class_pdb(oxygens, 4, h1_final, h2_final, 3, startn)
        oxygens_to_print = water_positions.oxy_processed()
        hydro1_to_print = water_positions.hydro_processed1()
        hydro2_to_print = water_positions.hydro_processed2()
        
        water = oxygens_to_print + hydro1_to_print + hydro2_to_print
        water_s = sorted(water)
        
#pieces of water:
        for j in water_s:
            new_line = '{:<0} {:>6} {:>4} {:>3} {:>5} {:>11} {:>7} {:>7} {:>5} {:>5} {:>11}'.format( \
                        j[0], j[1], j[2], j[3], j[4], j[5], j[6], j[7], j[8], j[9], j[10])
            print(new_line, file = waterout)
        insert_str = 'TER {} \n'.format(len(water_s)+len(atoms_tube)+1)
        waterout.write(insert_str)
            #waterout.write(new_line)
        
        conect = water_conections(oxygen_quantity, startn)
        for k in conect:
            print(*k, file = water_conect_out)
        water_conect_out.write('END')

        
#combining all together (merging) and getting the final file!
        #final_pdb_water_and_tube = header_part + atoms_tube + water_s + conect_tube + conect
    
        filenames = [header_out, tubeout, waterout, tube_conect_out, water_conect_out]
            
        for fname in filenames:
            fname.seek(0)
            #with open(fname) as infile:
            for line in fname:
                    final_output.write(line)
    
    '''This steps can be done before creating data_file, since in theory there are no *.txt files other than the ones
        created by this code '''
    for filename in glob.glob("./*.txt"):
        os.remove(filename)


class get_list_for_pdb_input_file(object):
    
    def __init__(self, file_get_center_tube):
        self._infile = file_get_center_tube

    def getting_Atoms(self):
        coordinatesxy = []
        for j, lines in enumerate(self._infile):
                
            if j > 8:
                elementspdb = lines.split()
                if elementspdb[0] =='ATOM': 
                    coordinatesxy.append(elementspdb)
        self._infile.seek(0)
        return coordinatesxy
    
    def getting_header_from_pdb_file(self):
        header_list = []
        for j, lines in enumerate(self._infile):
            if j <= 8: 
                elementspdb = lines.split()
                header_list.append(elementspdb)
        self._infile.seek(0)
        return header_list

    def getting_conect(self):
        conect_list = []
        for _j, lines in enumerate(self._infile):
                elementspdb = lines.split()
                if elementspdb[0] =='CONECT': 
                    conect_list.append(elementspdb)
        self._infile.seek(0)
        return conect_list
    
        


def gettin_center_of_tube(file_with_tube):
        coordinatesxy = []
        for j, lines in enumerate(file_with_tube):
            
            if j > 8:    #this works if the file only contains the nanotube, nothing else, not water nor whatever
                elementspdb = lines.split()
                if elementspdb[0] =='ATOM': 
                    coordsxy = [float(elementspdb[5]), float(elementspdb[6]), float(elementspdb[7])]
                    coordinatesxy.append(coordsxy)
             
        xyz1 = random.choice(coordinatesxy)
        xyz2 = random.choice(coordinatesxy)
        xyz3 = random.choice(coordinatesxy)

        while xyz1[0] == xyz2[0] or xyz1[0] == xyz3[0] or xyz2[0] == xyz3[0] or \
            xyz1[1] == xyz2[1] or xyz1[1] == xyz3[1] or xyz2[1] == xyz3[1]:
            
            xyz1 = random.choice(coordinatesxy)
            xyz2 = random.choice(coordinatesxy)
            xyz3 = random.choice(coordinatesxy)
        
        circle_points_3 = [xyz1, xyz2, xyz3]
        circle_points_3_s = sorted(circle_points_3, key = lambda list_of_3:list_of_3[0])
        #print(circle_points_3_s)
        
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
    while i < ((3*oxygen_quantity)+startn+1):
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
        self._atoms = 'ATOM'
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

def calculate_water_molecules_inside_nanotube(D, L, density_g_cc):
    water_mass = 18.01528
    avogadro = 6.022e23
    cm3_ti_ang3 = 1e24
    R = D/2.0
    v_ang3 = R*R*L*np.pi
    den_molecules_ang3 =    avogadro*density_g_cc/(water_mass*cm3_ti_ang3)
    water_molecules_to_fill = den_molecules_ang3*v_ang3
    return water_molecules_to_fill


def assign_positions_to_oxygens_list_pdb(Diameter, Length, number_of_atoms, bond_length, cylinder_x, cylinder_y, safety_clearance):
    L = Length - (4*bond_length) #safety clearance
    D = Diameter 
    list_positions_0 = []
    list_positions_final = []
    
    increase_by_fix = (1/(number_of_atoms/L)) 
    increase_by = 1.0
    r= (D/2.0)- safety_clearance #safety clearance
    
    while len(list_positions_0) <= number_of_atoms:
        #getting coordinates from (0,0,0)
        x = np.random.uniform(-1,1)*(r)
        y = np.random.uniform(-1,1)*(r)
        #z = (np.random.uniform(0,1)*L)+ (2*bond_length)
        r_test = np.sqrt(x*x + y*y)
        
        if r_test < r:
            increase_by += increase_by_fix
            #print(np.random.uniform(0,1))
            z = (np.random.uniform(0,1)*1)  + increase_by                        #+ (2*bond_length)
            coords_0 = [x,y,z]
            list_positions_0.append(coords_0)
            #L += 1.0
        #print(L)
    
    #changin coordinates, translation
    for xyz in list_positions_0:
        new_x = xyz[0] + cylinder_x
        new_y = xyz[1] + cylinder_y
        new_z = xyz[2]
        coords_final = [new_x, new_y, new_z]
        list_positions_final.append(coords_final)
    #print(list_positions_final)
    return list_positions_final

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

#     totalh1 = []
#     values = [0,0,0]
#     #sprint(values)
#     for _numbers_atoms in range((1*number_of_atoms)+1):
#         u = np.random.uniform()
#         v = np.random.uniform()
#         theta = 2*np.pi*u 
#         phi = mt.acos(2*v-1)
#         xh1 = values[0] + (bond_length*np.sin(phi)*np.cos(theta))
#         yh1 = values[1] + (bond_length*np.sin(phi)*np.sin(theta))
#         zh1 = values[2] + (bond_length*np.cos(phi))
#         h1 = [xh1, yh1, zh1]
#         totalh1.append(h1)
#     return totalh1


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
    
    
#     totalh2 = []
#     for lines in h1_list:
#         xh1 = lines[0]
#         yh1 = lines[1]
#         zh1 = lines[2]
#         xyz = np.random.randint(1, 4)
#         #xyz = 1
#         if xyz == 3: #rotates in z
#             xh2 = xh1*np.cos(angle) - yh1*np.sin(angle)
#             yh2 = xh1*np.sin(angle) + yh1*np.cos(angle)
#             zh2 = zh1
#          
#         elif xyz == 2: #rotates in y
#             xh2 = xh1*np.cos(angle) + zh1*np.sin(angle)
#             yh2 = yh1
#             zh2 = -xh1*np.sin(angle) + zh1*np.cos(angle) 
#             
#         if xyz == 1: #rotates in x
#             xh2 = xh1
#             yh2 = (yh1*np.cos(angle) - zh1*np.sin(angle))
#             zh2 = (yh1*np.sin(angle) + zh1*np.cos(angle))
#         
#         
#         h2 = [xh2, yh2, zh2]
#         totalh2.append(h2)
#         
#     return totalh2

################   END    
if __name__ == "__main__": main()
