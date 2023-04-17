import ase
import sys
import pandas as pd
import itertools
import os, getpass
from pymatgen.core.composition import Composition
from mpds_client import MPDSDataRetrieval, MPDSDataTypes, APIError

def is_disorder(formula, occ, logfile):
    if occ is None:
        return False
    elif not all([i.is_integer() for i in occ]):
        with open(logfile, 'a') as f:
            print(formula, file=f)
        return True

def symbols(formula):
    try:
        CC = Composition(formula)
    except:
        formula = ''.join(filter(str.isalpha, formula))
        CC = Composition(formula)
    c = [i for i in CC.as_dict().keys()]
    return set(c)

def wrong_contents(elements, formula):
    """ fit description"""
    els_mpds = symbols(formula)
    #els_mpds = set(formula)
    if len(els_mpds) > len(elements) or \
            len(els_mpds.difference(elements)):
        print('wrong content:', formula)
        return True
    else:
        return False

def combinations(elements,n):
    return ['-'.join(list(i)) for i in set(itertools.combinations(elements, n))]

def pullMPDS(client, qe, clas, dirname, logfile):
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

def get_references(elements, dirname, logfile='mpds_reference_log'):
    os.environ['MPDS_KEY'] = 'ggBYYU0tszpYMTqLahr604WPM3Ao8o5lK3XTCV46FjyR0j2y'
    client = MPDSDataRetrieval(dtype=MPDSDataTypes.PEER_REVIEWED)

    _elements = set(elements.split('-'))
    classes = {1:"unary",2:"binary",3:"ternary",4:"quaternary", 5:"quinary"}

    for n, clas in classes.items():
        query_elements = combinations(_elements, n) 
        for qe in query_elements:
            print(f"{clas}", qe)
            try:
                answer = pullMPDS(client, qe, clas, dirname, logfile)
            except:
                continue
            for item in answer:
                #if is_disorder(item[0], item[-1], f'{dirname}/{logfile}'):
                #    print('disordered:', item[0])
                #    continue
          
                 crystal = MPDSDataRetrieval.compile_crystal(item, 'ase')
          
                 if not crystal:
                     print(f'Cannot build crystal')
                     continue
                 else:
                     print('Crystal seems right', f'{item[2]}_{item[0]}')

                 rename_crystal(crystal, item[0], dirname)


if __name__=="__main__":
    elements = ['Mg', 'Al', 'Cl', 'Sn']
    elements = '-'.join(elements)
    dirname = elements
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    get_references(elements, dirname=dirname, logfile=f'mpds_references_{elements}_log')
