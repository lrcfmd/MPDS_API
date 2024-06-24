import ase
from ase import neighborlist
import sys
import pandas as pd
import os, getpass
from mpds_client import MPDSDataRetrieval, MPDSDataTypes, APIError

def calculate_neighbors(dic, mol):
    cut = neighborlist.natural_cutoffs(mol)
    nl = neighborlist.NeighborList(cut, self_interaction=False, bothways=True)
    nl.update(mol)

    nn = {i.symbol : set() for i in mol if i.symbol in dic}

    for i, n in enumerate(mol):
        ind, _ = nl.get_neighbors(i)
        if n.symbol in dic:
            nn[n.symbol].add(len(ind))

    for n, l in nn.items():
        if n in dic: 
            dic[n].extend(list(l))
        else:
            dic[n] = list(l)

    return dic

def calculate_lengths(ase_obj, elA, elB, limit=13):
    lengths = []
    for n, atom in enumerate(ase_obj):
        if atom.symbol == elA:
            for m, neighbor in enumerate(ase_obj):
                if elA == elB and m == n: continue
                if neighbor.symbol == elB:
                    try:
                        dist = round(ase_obj.get_distance(n, m), 2) # NB occurrence <-> rounding
                    except: continue
                    if dist < limit:
                        lengths.append(dist)
    return list(set(lengths))

os.environ['MPDS_KEY'] = 'ggBYYU0tszpYMTqLahr604WPM3Ao8o5lK3XTCV46FjyR0j2y'
datatypes = [x for x in dir(MPDSDataTypes) if not x.startswith('__')]

client = MPDSDataRetrieval(dtype=MPDSDataTypes.PEER_REVIEWED)

sets = ["Y-Sr-O", "Ti-Sr-O", "Y-Ti-O", "Y-O", "Sr-O","Ti-O" ]
dfrm = pd.DataFrame(columns=['O-O'])
lengths = []

for s in sets:
    answer = client.get_data(
        {"elements": s, "props": "atomic structure", "classes": "non-disordered"},
        fields={'S':[
                'phase_id',
                'entry',
                'chemical_formula',
                'cell_abc',
                'sg_n',
                'basis_noneq',
                'els_noneq'
                ]}
    )
    for item in answer:
        crystal = MPDSDataRetrieval.compile_crystal(item, 'ase')
        if not crystal: continue

        lengths.extend(calculate_lengths(crystal, 'O', 'O'))

dfrm['O-O'] = sorted(lengths)
dfrm.to_csv(f'SrTiOY_bondlenghts_perN/sortedhO-O.csv', index=None)
dfrm['occurrence '+ 'O-O' ] = dfrm.groupby('O-O')['O-O'].transform('count')
dfrm.drop_duplicates('O-O', inplace=True)

dfrm.to_csv(f'SrTiOY_bondlenghts_perN/bondLengthO-O.csv', index=None)
