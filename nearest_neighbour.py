# A script to extract info of nearest neighbours and connectivity of atoms 
# from the structure

import sys
from ase.io import read
import matplotlib.pyplot as plt
from ase import neighborlist
import numpy as np
from sklearn.metrics.pairwise import euclidean_distances
#from ase.geometry.analysis import Analysis
#from scipy import sparse

def name(i, j, mol):
    return ' '.join(list(np.sort([mol[i].symbol, mol[j].symbol])))

def update_distances_dic(mol, i, j, d):
    """ dic has to be instantiated with all items"""
    distance_dic[name(i, j, mol)] += [d]

def update_neighbors_dic(symbol, nn):
    """ dic has to be instantiated with all items"""
    neighbors_dic[symbol] += [d] 

mol = read('Br12Ca1Li2Zr2.cif')

cut = neighborlist.natural_cutoffs(mol)
#cut ={'Br':1.2, 'Ca':1.76, 'Li':1.4, 'Zr': 1.75}

# neighbours list
nl = neighborlist.NeighborList(cut, self_interaction=False, bothways=True)
nl.update(mol)

# nearest neighbours

#nnei = {i:[] for i in mol.get_chemical_symbols()}

for i, n in enumerate(mol):
    ind, displacements = nl.get_neighbors(i)
    #nnei[n.symbol].append(len(ind))
    update_neighbors_dic(n.symbol, len(ind))

    for j,d in zip(ind,displacements):
        nei_position = mol.positions[j] + d @ mol.get_cell()
        d = euclidean_distances([mol.positions[i]], [nei_position])[0][0]) 
        update_distances_dic(mol, i, j, d)


nn = {i:set(n) for i, n in nnei.items()}

for i, n in nn.items():
    print(f"For Element {i} # nearest neighbours: {n}")


sys.exit()
# NN distances:

def make_dic(i, j, d):
    """ dic has to be instantiated with all items"""
    dic[name(i,j)] += [d]


i,j,d = neighborlist.neighbor_list('ijd', mol, cut)

dic = { name(i,j) : [] for i, j in zip(i,j)}

for i, j, d in zip(i,j,d):
    make_dic(i,j,d)

for k,v in dic.items():
    print(k, np.average(v))


# Connectivity
#matrix = neighborlist.get_connectivity_matrix(nl.nl, sparse=False)
#connected_distance = neighborlist.get_distance_matrix(matrix)
#print(connected_distance)

#n_components, component_list = sparse.csgraph.connected_components(matrix)
#print(n_components, component_list)
