mmport pandas as pd
import ase.io.vasp
import ase.io
import sys

filename=sys.argv[1]
atoms=ase.io.read(filename)

cell = atoms.get_cell()
symbols = atoms.get_chemical_symbols()
positions = [tuple(i) for i in atoms.get_positions()]
formula = atoms.get_chemical_formula()

df = pd.DataFrame({'atoms':symbols, 'positions':positions})
df = df.sort_values(['atoms'])

batoms = ase.Atoms(formula, cell=cell)
batoms.set_chemical_symbols(list(df.atoms.values))
batoms.set_positions(list(df.positions.values))
print(formula)
batoms.write(f'{formula}_POSCAR', format='vasp')


#ase.io.write('input.cif', atoms)
#ase.io.vasp.write_vasp('TEST_POSCAR', atoms, direct=True)
