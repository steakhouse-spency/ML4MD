
from helperfun import *




'''
    parse pdb with only nanotube atoms/conects

'''
if debug:
    print("label: "+ label)

# open pdb file
empty_tube_file = open("../empty_tube/{}.pdb".format(label), 'r')
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
water_conect = water_conections(len(oxygens), startn)


centerAtoms(tube_atoms, tube_center)
centerAtoms(water_atoms, tube_center)

'''
    output pdb of water-filled nanotube

'''
if debug:
    print("new tube center: ", get_tube_center(tube_atoms))
    print("tube atoms: %d" % len(tube_atoms))
    print("water atoms: %d" % len(water_atoms))
    print("total atoms: %d" % (len(water_atoms) + len(tube_atoms)))

# final_out = open("../filled_tube/{}.pdb".format(label), 'w')
# rebo/norebo in filename
final_out = open("../filled_tube/{}-{}.pdb".format(label,sim), 'w')



final_out.write(header)

for atom in tube_atoms:
	line = formatPdbRow(atom, True)
	final_out.write(line)
for atom in water_atoms:
    line = formatPdbRow(atom, False)
    final_out.write(line)
for conect in tube_conect+water_conect:
    final_out.write(conect)

final_out.write("END\n")
final_out.close()




'''
    output datafile of water-filled nanotube

'''

bonds = []
angles = []

# with open('data_file/{}.data'.format(label), 'w') as 
# datafile = open('../data_file/{}.data'.format(label), 'w')
# rebo/norebo
datafile = open('../data_file/{}-{}.data'.format(label,sim), 'w')

# format bonds,angles for water only
for conect in water_conect:
    conect = conect.split()
    if(len(conect) == 4):
        formatWaterBond(bonds, conect)
        formatWaterAngle(angles, conect)

dataFileHeader(datafile, len(tube_atoms)+len(water_atoms), len(bonds), len(angles))
dataFileMass(datafile)
dataFileAtoms(datafile, tube_atoms+water_atoms)
dataFileBonds(datafile, bonds)
dataFileAngles(datafile, angles)