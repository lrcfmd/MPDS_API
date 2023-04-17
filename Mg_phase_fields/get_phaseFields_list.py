import ase
import sys
import pandas as pd
import itertools
import os, getpass
from pymatgen.core.composition import Composition
from mpds_client import MPDSDataRetrieval, MPDSDataTypes, APIError

def combinations(elements,n):
    return ['-'.join(list(i)) for i in set(itertools.combinations(elements, n))]

def pullMPDS(client, qe, clas):
    return client.get_data({"elements": qe, "props": "atomic structure", "classes": clas},
                           fields={'S':[
                                   'phase_id',
                                   'entry',
                                   'chemical_formula',
                                   #'chemical_elements',
                                   'cell_abc',
                                   'sg_n',
                                   'basis_noneq',
                                   'els_noneq'
                                   #'occs_noneq'
                                   ]}
                          )


def rename_crystal(atoms, ids, dirname):
    """ reorder atoms by alphabet; rename crystal by formula """
    cell = atoms.get_cell()
    symbols = atoms.get_chemical_symbols()
    positions = [tuple(i) for i in atoms.get_positions()]
    formula = atoms.get_chemical_formula()
    
    df = pd.DataFrame({'atoms':symbols, 'positions':positions})
    df = df.sort_values(['atoms'])
    
    batoms = ase.Atoms(formula, cell=cell)
    batoms.set_chemical_symbols(list(df.atoms.values))
    batoms.set_positions(list(df.positions.values))
    batoms.write(f'{dirname}/{formula}_{ids}_POSCAR', format='vasp')

def get_references(ids, elements, dirname, logfile='mpds_reference_log'):
    os.environ['MPDS_KEY'] = 'ggBYYU0tszpYMTqLahr604WPM3Ao8o5lK3XTCV46FjyR0j2y'
    client = MPDSDataRetrieval(dtype=MPDSDataTypes.PEER_REVIEWED)

    _elements = set(elements.split('-'))
    classes = {1:"unary",2:"binary",3:"ternary",4:"quaternary", 5:"quinary"}


    for n, clas in classes.items():
        query_elements = combinations(_elements, n) 
        for qe in query_elements:
            print(f"{clas}", qe)
            try:
                answer = pullMPDS(client, qe, clas)
            except:
                continue
            for item in answer:
                 crystal = MPDSDataRetrieval.compile_crystal(item, 'ase')
          
                 if not crystal:
                     print(f'Cannot build crystal')
                     continue
                 elif f'{item[2]}_{item[0]}' not in ids:
                     ids.append(f'{item[2]}_{item[0]}')
                     print('Crystal seems right', f'{item[2]}_{item[0]}')

                     rename_crystal(crystal, item[0], dirname)


if __name__=="__main__":
    phases = pd.read_csv('MgMMA_top20.csv')
    phases = [i.split() for i in phases['phases']]
    ids = []
    for elements in phases:
        elements = '-'.join(elements)
        dirname = elements
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        get_references(ids, elements, dirname=dirname, logfile=f'mpds_references_{elements}_log')
    pd.DataFrame({'ids':ids}).to_csv('Downloaded_ids.csv', index=False)
