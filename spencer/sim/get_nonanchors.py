
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

def get_nonanchors(tube_z):

	# number of tube atoms
	tube_len = len(tube_z)

	# get sorted indicies
	sorted_i = [i[0] for i in sorted(enumerate(tube_z), key=lambda x:x[1])]
	
	for i in range(tube_len-1, 0, -1):
		print(
		
	
	# where to place the ring
	rings=[0.25, 0.5, 0.75, 1.0]

	# atoms per ring
	ring_len = 10
	half = int(ring_len/2)

	# for every ring
	start = ring_len
	rng = [None, None]
	non_anchor_i = []
	anchors = [*range(0, ring_len)] + [*range(tube_len-ring_len, tube_len)]
	for ring in rings:	
		mid = int(tube_len*ring)
		# get range of atoms in ring
		rng[0] = mid - (half - 1)
		rng[1] = mid + half + 1

		anchors += [*range(rng[0], rng[1])]
	print("anchors:")
	print(anchors)
	exit(1)

'''	for ring in rings:

		# get id of atom in middle of ring
		mid = int(tube_len*ring)
		
		# get range of atoms in ring
		# using 10 atoms/ring 
		rng[0] = mid - (half - 1)
		rng[1] = mid + half

		anchor += [*range(rng[0], rng[1]+1)]
		


		# append non-anchor atom ids
		if ring != 1.0:
			# from 'start' to the first atom in this ring
			end = rng[0]
		else:
			# from end of last ring to end of tube
			# -10 for last ring at end of tube
			#echo "last iter"
			end = tube_len - ring_len
			print(end, tube

		# append new range of atom ids
		non_anchor_i += [*range(start, end)]

		# update start to the last atom in this ring + 1
		start = rng[1] + 1
'''
	# convert indicies to atom ids
	non_anchor_id = []
	for i, anchor_i in enumerate(non_anchor_i):
		non_anchor_id.append(sorted_i[anchor_i] + 1)
		# sorted_i[i] += 1

	# sort and return ids
	non_anchor_id.sort()

	return non_anchor_id

tube_z, water_id = getAtomData()
tube_id = get_nonanchors(tube_z)

print(len(tube_id))
print(len(water_id))
for i in tube_id:
	print(i)


