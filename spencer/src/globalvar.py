from nanotube_radius import getDiam

# global variables
debug = True


# nano tube structure parameters
n = 7
m = 14

## this is is nm not ang... okay?
D = getDiam(n, m)
if debug:
	print("tube diameter: %f nm" % D)

#nanotube length (ang)
L = 250

# nano tube material
material = "C"


density = 0.471



# pdb file path
# only for empty nano tubes
# pdb must only have: 1 line header, nano tube atoms of type 'residue', and CONECT for nano tube atoms
label = "{}_{}_{}_{}".format(material, int(L/10), n, m)




# DONT CHANGE, same values for every simulation
#############
box = 150
den_mol = 29.91577548987 #A^3/molecule
side_cube = den_mol**(1/3)

#density:
rho = 0.471 #value obtained from previous simulation. This is a pseudo density
cutoff = 8.0

sigma_O = 3.1668
eps_O = 0.21084

#water parameters:
bond_length = 0.9572
angle = 104.52
safety_clearance = 1.0


# residue within pdb file for atoms (e.g. HETATM, ATOM, etc.)
residue = "XXX"
record_types = ["HETATM", "ATOM"]


elem_types = ['C', 'H', 'O']

elem_mass = {'Si': 28.085,
			 'C': 12.011, 
			 'H': 1.008, 
			 'O': 15.9994}

elem_charge = {	'C': 0.0,
				'H': 0.5242, 
			    'O': -1.0484}
