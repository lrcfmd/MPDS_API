import pandas as pd
from pymatgen.core.composition import Composition as C

df = pd.read_csv('mpds_scon.csv')

formulas = df.values[:,2]
uniq = []
for f in formulas:
    try:
        f = C(f.split()[0]).as_dict()
        f = ' '.join(list(f.keys()))
    except:
        pass
    if f not in uniq:
         uniq.append(f)
         print(f)

print(len(formulas), len(uniq)) 
