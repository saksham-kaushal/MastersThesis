# import numpy as np

import core

# =============================== User-defined functions ======================

# =============================== Main program ======================

if __name__ == '__main__':

	assembly_list 			= ['GM-Early','Organic','GM-Late']
	# assembly_list 				= ['Organic']

	for assembly in assembly_list:
		# ------- Get HDF5 file handles list for each type of assembly and store in a dictionary.

		files  = {core.get_index(fname.filename.split('/')[-1]):core.get_files(core.get_directory(str(assembly)),fname.filename.split('/')[-1]) for fname in core.get_files(core.get_directory(str(assembly)))}
		print(files)
		print(core.get_index_redshift_dict())