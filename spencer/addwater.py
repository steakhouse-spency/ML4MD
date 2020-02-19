
from helperfun import *




'''
    parse pdb with only nanotube atoms/conects

'''

# open pdb file
empty_tube_file = open(input_file, 'r')
# seperate empty tube pdb file by header, atoms, conect
header, tube_atoms, tube_conect = seperate_empty_tube(empty_tube_file)
# close it
empty_tube_file.close()


'''
    compute pdb atoms/conects for water molecules
    
'''
# get x,y coordinants representing the center of the Nano Tube  
tube_center = get_tube_center(tube_atoms)
startn = len(tube_atoms)+1
oxygen_quantity = calculate_water_molecules_inside_nanotube(rho)
oxygens = assign_positions_to_oxygens_list_pdb_new(oxygen_quantity, tube_center, startn)
hydro1 = adding_H1(len(oxygens), angle, bond_length)
hydro2 = adding_H2(hydro1, angle, bond_length)
h1_final = final_position_for_H(oxygens, hydro1) #move the water atoms from the origin to the random positions generated in oxygens
h2_final = final_position_for_H(oxygens, hydro2)
water_positions = water_class_pdb(oxygens, 4, h1_final, h2_final, 3, startn)
oxygens_to_print = water_positions.oxy_processed()
hydro1_to_print = water_positions.hydro_processed1()
hydro2_to_print = water_positions.hydro_processed2()
water_atoms = oxygens_to_print + hydro1_to_print + hydro2_to_print
water_atoms = sorted(water_atoms)
water_conect = water_conections(oxygen_quantity, startn)



'''
    output pdb of water-filled nanotube

'''

final_list = tube_atoms + water_atoms + tube_conect + water_conect
final_out = open("output.pdb", 'w')
final_out.write(header)
# convert every row to one string and write to final pdb
for row in final_list:
    line = ""    
    for i, col in enumerate(row):
        if i >= 6 and i <= 8:
            # print(col)
            # print(round(col,3))
            rounded =  str(round(float(col),3))
            line += rounded + "\t"
        else:
            line += str(col) + "\t"

    line = line[:-1] + "\n"
    final_out.write(line)
final_out.write("END\n")
final_out.close()




'''

    Things to do:

        - normalize coordinates?
        - check assign_positions_to_oxygens_list_pdb_new


    Extra:

        - make code nicer via some oop?


'''