##Script for pulling reference structures from MPDS

Script can be run iteratively for a list of phase fields.
In this case, the library of downloaded structures is created to avoid duplicates.
The results are saved in the respective folders and summarized in 'Downloaded_ids.csv' file.

1. Prepare a list of phase fields (e.g., ranked list with VAE, 'MgMMA_top20.csv') 
2. Modify get_Phasefields_list.py file:
a. line X to include your MPDS API user key
b.  line 73 to address the list of phase fields
3. run ' python get_Phasefields_list.py' 

## TODO
# Flag disordered structures
