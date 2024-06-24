import ase
from ase import neighborlist
import sys
import pandas as pd
import os, getpass
from mpds_client import MPDSDataRetrieval, MPDSDataTypes, APIError

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

dfrm = pd.DataFrame(columns=['O-O'])
lengths = []

for s in ['Li-O-Mg-P']: #, 'Mg-O', 'P-O','Al-O', 'Li-O-Al', 'Li-O-P', 'Li-O-Mg', 'Mg-Al-O', 'Mg-P-O']: #, 'Al-P-O']:
    answer = client.get_data(
        {"elements": s, "props": "atomic structure"}, #"classes": "non-disordered"},
        fields={'S':[
#                'phase_id',
#                'entry',
                'chemical_formula',
#                'cell_abc',
#                'sg_n',
#                'basis_noneq',
#                'els_noneq',
                'condition'
                ]}
    )
    for item in answer:
        if item[-1] and item[-1][0] and (290 < item[-1][0] < 300):
            print(item[-1][0])

#       crystal = MPDSDataRetrieval.compile_crystal(item, 'ase')
#       if not crystal: continue
#       lengths.extend(calculate_lengths(crystal, 'O', 'O'))
sys.exit()
dfrm['O-O'] = sorted(lengths)
dfrm.to_csv(f'sortedhO-O.csv', index=None)
dfrm['occurrence '+ 'O-O' ] = dfrm.groupby('O-O')['O-O'].transform('count')
dfrm.drop_duplicates('O-O', inplace=True)

dfrm.to_csv(f'bondLengthO-O.csv', index=None)
