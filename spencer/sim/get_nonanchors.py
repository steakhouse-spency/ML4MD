
import os
import sys
import subprocess

label=sys.argv[1]

material = "C"

def getAtomData():

	fpdb = open("%s.pdb" % (label), "r")
	tube_z = []
	water_id = [] 
	for line in fpdb:
		cols = line.split()
		if cols[0] == "ATOM":
			if cols[2] == material:
				tube_z.append(float(cols[7]))
			elif cols[2] == "MOL":
				water_id.append(cols[1])
			else:
				print("error")
				print(">%s<" % (cols[2]))
				exit(1)

	return tube_z, water_id

def Diff(li1, li2): 
    return (list(set(li1) - set(li2))) 

def index_to_id(i_list, i_map):
	ids = []
	for i in i_list:
		ids.append(i_map[i]+1)
	return ids

def get_nonanchors(tube_z):

	# number of tube atoms
	tube_len = len(tube_z)

	# get sorted indicies
	sorted_i = [i[0] for i in sorted(enumerate(tube_z), key=lambda x:x[1])]				
	
	# where to place the ring
	rings=[0.25, 0.5, 0.75]

	# atoms per ring
	ring_len = 10
	half = int(ring_len/2)

	# for every ring
	# start = ring_len
	rng = [None, None]
	# non_anchor_i = []
	anchors = [*range(0, ring_len)] + [*range(tube_len-ring_len, tube_len)]
	
	for ring in rings:	
		mid = int(tube_len*ring)
		# get range of atoms in ring
		rng[0] = mid - (half - 1)
		rng[1] = mid + half + 1

		anchors += [*range(rng[0], rng[1])]
	

	anchors = index_to_id(anchors, sorted_i)

	all_atoms = [*range(tube_len)]
	all_atoms = index_to_id(all_atoms, sorted_i)

	non_anchor = Diff(all_atoms, anchors)

	return non_anchor

tube_z, water_id = getAtomData()
non_anchor = get_nonanchors(tube_z)

print(len(tube_z))
print(len(water_id))
for i in non_anchor:
	# increment by 1
	print(i)


