# global variables


# DONT CHANGE, same values for every simulation
#############
box = 150
den_mol = 29.91577548987 #A^3/molecule
side_cube = den_mol**(1/3)

#density:
rho = 0.471 #value obtained from previous simulation. This is a pseudo density
cutoff = 8.0

#nanotube dimensions
L = 82.817
sigma_O = 3.1668
eps_O = 0.21084

#water parameters:
bond_length = 0.9572
angle = 104.52
safety_clearance = 1.0


# residue within pdb file for atoms (e.g. HETATM, ATOM, etc.)
residue = "XXX"
record_types = ["HETATM", "ATOM"]

# nano tube structure parameters
m = 15
n = 7
D = 15.2
density = 0.471
material = "C"

# pdb file path
# only for empty nano tubes
# pdb must only have: 1 line header, nano tube atoms of type 'residue', and CONECT for nano tube atoms
label = "{}_{}".format(material, D)

elem_types = ['C', 'H', 'O']

elem_mass = {'Si': 28.085,
			 'C': 12.011, 
			 'H': 1.008, 
			 'O': 15.994}

elem_charge = {'Si': 0.6,
			   'C': -0.6, 
			   'H': 0.5897, 
			   'O': -1.1794}
