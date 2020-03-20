import collections 
from collections import Counter
import pandas as pd
import numpy as np

def main():
    myfile = "thermalizedsicnt-end-298-0-min.SiCNT"
    line_where_Bonds_start = 5298  #end of the atoms section - 1 (as read in eclipse editor)
    with open('new_data_file.data', 'w') as outfile, open(myfile, 'r') as infile, \
         open('data-file-sicnt-water-0.6.data', 'r') as firstdatafile:
        atom_list_new = positions_list(infile)
        
        #lines = f.readlines()
        #lines = [l for l in lines if "ROW" in l]
        #outfile.writelines(lines)
        
        #get the header of the first data file
        for j, linei in enumerate(firstdatafile):
            if j < 25: #where header end, this one is not gonna change
                outfile.write(linei)
        firstdatafile.seek(0)
        
        
        #print the new coordinates:
        for i in atom_list_new:
            print(*i, file = outfile)
        
        
        #get the bonds and the angles
        for j, linei in enumerate(firstdatafile):
            if j > line_where_Bonds_start:
                outfile.write(linei)
        
        
        
        
        
def get_header_from_first_data_file():
        pass
    
def positions_list(file_tx): #return ordered atoms 
    dict_atoms = {}
    coords = []
    for skip in range(9):
        next(file_tx)
    for rows in file_tx:
            values = rows.split()
            keys = int(values[0])
            atom_type = int(values[2])
            if atom_type == 1: #Si
                q = 0.6
                mol_id = 1
            elif atom_type == 2: #C
                q = -0.6
                mol_id = 1
            
            elif atom_type == 3: #H
                q = 0.5897
                mol_id = 2
            
            elif atom_type == 4: #O
                q = -1.1794
                mol_id = 2
                
            coords = [mol_id, int(values[2]), q, float(values[3]), float(values[4]), float(values[5])]
            
            
            dict_atoms[keys] = coords   
    positions = sorted(dict_atoms.items())    
    file_tx.seek(0)
    
    #convert to a list:
    final_list = []
    for i in positions:
        ncoords = [i[0], i[1][0], i[1][1], i[1][2], i[1][3], i[1][4], i[1][5]]
        final_list.append(ncoords)
        
        
    return final_list
        
if __name__ == "__main__": main()
