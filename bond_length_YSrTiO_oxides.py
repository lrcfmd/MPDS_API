import os
import pandas as pd
from ase import neighborlist
from mpds_client import MPDSDataRetrieval, MPDSDataTypes, APIError

def calculate_neighbors(atom_dict, molecule):
    """
    Calculate the number of neighbors for each atom type in the molecule.

    Parameters:
    atom_dict (dict): Dictionary to store the neighbor counts.
    molecule (ase.Atoms): Molecule object.

    Returns:
    dict: Updated dictionary with neighbor counts.
    """
    cutoffs = neighborlist.natural_cutoffs(molecule)
    nl = neighborlist.NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(molecule)

    neighbor_counts = {atom.symbol: set() for atom in molecule if atom.symbol in atom_dict}

    for i, atom in enumerate(molecule):
        indices, _ = nl.get_neighbors(i)
        if atom.symbol in atom_dict:
            neighbor_counts[atom.symbol].add(len(indices))

    for symbol, counts in neighbor_counts.items():
        if symbol in atom_dict:
            atom_dict[symbol].extend(list(counts))
        else:
            atom_dict[symbol] = list(counts)

    return atom_dict

def calculate_lengths(ase_obj, element_a, element_b, limit=13):
    """
    Calculate the bond lengths between two elements within a specified limit.

    Parameters:
    ase_obj (ase.Atoms): ASE atoms object.
    element_a (str): Symbol of the first element.
    element_b (str): Symbol of the second element.
    limit (float): Maximum distance to consider for bond lengths.

    Returns:
    list: Unique bond lengths between element_a and element_b within the limit.
    """
    lengths = []
    for i, atom_a in enumerate(ase_obj):
        if atom_a.symbol == element_a:
            for j, atom_b in enumerate(ase_obj):
                if element_a == element_b and j == i:
                    continue
                if atom_b.symbol == element_b:
                    try:
                        distance = round(ase_obj.get_distance(i, j), 2)
                    except Exception:
                        continue
                    if distance < limit:
                        lengths.append(distance)
    return list(set(lengths))

# Set MPDS API key
os.environ['MPDS_KEY'] = 'your_mpds_key'  # INSERT_KEY

# Initialize MPDS client
client = MPDSDataRetrieval(dtype=MPDSDataTypes.PEER_REVIEWED)

# Define sets of elements to search for
element_sets = ["Y-Sr-O", "Ti-Sr-O", "Y-Ti-O", "Y-O", "Sr-O", "Ti-O"]
dataframe = pd.DataFrame(columns=['O-O'])
all_lengths = []

# Retrieve data and calculate bond lengths
for elements in element_sets:
    try:
        response = client.get_data(
            {"elements": elements, "props": "atomic structure", "classes": "non-disordered"},
            fields={'S': [
                'phase_id',
                'entry',
                'chemical_formula',
                'cell_abc',
                'sg_n',
                'basis_noneq',
                'els_noneq'
            ]}
        )
    except APIError as e:
        print(f"API Error: {e}")
        continue

    for item in response:
        crystal = MPDSDataRetrieval.compile_crystal(item, 'ase')
        if not crystal:
            continue

        all_lengths.extend(calculate_lengths(crystal, 'O', 'O'))

# Process and save results
dataframe['O-O'] = sorted(all_lengths)
dataframe.to_csv('SrTiOY_bondlengths_perN/sorted_O-O.csv', index=False)
dataframe['occurrence O-O'] = dataframe.groupby('O-O')['O-O'].transform('count')
dataframe.drop_duplicates('O-O', inplace=True)
dataframe.to_csv('SrTiOY_bondlengths_perN/bondLength_O-O.csv', index=False)
