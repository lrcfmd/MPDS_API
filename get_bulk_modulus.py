import sys
import pandas as pd
import os, getpass
from mpds_client import MPDSDataRetrieval, MPDSDataTypes, APIError

os.environ['MPDS_KEY'] = 'your_mpds_key'  # INSERT_KEY
datatypes = [x for x in dir(MPDSDataTypes) if not x.startswith('__')]
example_props = [
#'electrical resistance',
#'Seebeck coefficient',
#'heat capacity at constant pressure',
'bulk modulus',
#'isothermal bulk modulus',
#'values of electronic band gap',
#'temperature fro congruent metling',
#'Debye temperature',
#'linear thermal expansion coefficient'
]

example_fields = {
    'P':[ # *P*hysical property entries
        'sample.material.entry',
        'sample.material.phase',
        'sample.material.chemical_elements',
        'sample.material.chemical_formula'
    ],
    'S':[ # Crystalline *S*tructure entries
        'entry'
        'phase',
        'chemical_elements',
        'chemical_formula'
    ],
    'C':[ # Phase diagrams, i.e. *C*onstitution entries
        'entry',
        lambda: 'MANY-PHASE', # constants are given like this (on purpose)
        'chemical_elements',
        lambda: 'MANY-FORMULAE'
    ]
    # NB. P-S-C are interconnected by means of the distinct phases
}

#myfield = {'P':['sample.material.chemical_formula', 'sample.material.chemical_elements'], 'S':['chemical_elements'], 'C':['chemical_elements']}
myfield = {'P':['sample.material.formula']}

client = MPDSDataRetrieval(dtype=MPDSDataTypes.MACHINE_LEARNING)
#client = MPDSDataRetrieval() #dtype=MPDSDataTypes.PEER_REVIEWED)
df = client.get_dataframe({"props": "isothermal bulk modulus"})
#df = client.get_dataframe({"classes": "superconductor"})
df.head()
df.to_csv('mpds_bulk_modulus.csv')
sys.exit(0)

#print(client.get_data({"elements": "O", "classes": "binary", "sgs": "I4/mmm"}))
#for entry in client.get_data({"classes": "superconductor"}):
#    print(entry, file=open('results_peer_reviewed_superconductor', 'a'))
#sys.exit(0)

scs = []

for card in client.get_dataframe({"props": "isothermal bulk modulus"}, fields=myfield):
     elements = ' '.join(card[0])
     if elements not in scs:
          scs.append(elements)
          print(f'elements,', file=open('mpds_peer_reviewed_superconductor.csv', 'a'))

print(f"Found {len(scs)} results")
