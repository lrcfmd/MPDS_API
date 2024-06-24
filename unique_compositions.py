import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re 
from pymatgen.core.composition import Composition as C

def findatom(f):
   f = re.split("[0-9\.\]\[\-\+\=\-\!\s]+", f)
   for s in ['rt', 'ht', 'hp', 'lt', 'hex', 'hom', 'orth', 'tet', 'tm', 'x', 'hyp', 'xlt', 'z']:
       if s in f: f.remove(s)
   return ''.join(f)

def parse_phases(dbfile, n, nt):

    print(f"Parsing {dbfile}...")
    df = pd.read_csv(dbfile)

    formulas = [findatom(f) for f in df.values[:,n]]
    temperes = [t for t in df.values[:, nt]]
    unique_compositions = []
    unique_t = []
    phases = []
    cant =[]

    dicT = {}

    for f, t in zip(formulas, temperes):
        try:
           pf = C(f).reduced_formula
        except:
            print (f"Cannot parse {f}")
            if f not in cant:
                cant.append(f)
            continue

        if pf not in list(dicT.keys()) or dicT[pf] < t:
             dicT[pf] = t

    #dicT = {f: t for f,t in zip(unique_compositions, unique_t)}
    dics = {f: t for f,t in sorted(dicT.items(), key=lambda x: x[0])}

    print('Unique compositions:', len(unique_compositions))
    print('Unique compositions couldnt parse:', len(cant))
    return set(unique_compositions), dics

#ferro, ferroT = parse_phases('mpds_ferromagnet.csv', 2, -1)
#aferro, aferroT = parse_phases('mpds_antiferromagnet.csv', 2, -1)

bulk, bulkM = parse_phases('mpds_bulk_modulus_peer.csv', 2, -1 )

plt.hist(bulkM.values(), bins='fd', label='Compositions with reported peer-reviewed bulk modulus')
#plt.hist(aferroT.values(), bins=40, label='Antiferromagnet', alpha=0.7)
plt.xlim(-10, 700)
plt.legend()
#plt.xlabel(r"T$_{magnetic transition}$, K", fontsize=14)
plt.xlabel(r"Isothermal bulk modulus, GPa", fontsize=14)
plt.ylabel("Counts", fontsize=14)
plt.show()
