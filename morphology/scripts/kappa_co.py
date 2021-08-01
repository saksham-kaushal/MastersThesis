import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import seaborn as sns

import core

# =============================== User-defined functions ======================

def get_angular_momentum(rxv, mass):
	return np.transpose(np.multiply(mass, np.transpose(rxv)))

# =============================== Main program ======================

if __name__ == '__main__':

	assembly_list 			= ['GM-Early','Organic','GM-Late']
	# assembly_list 				= ['Organic']
	kappa_co = dict()

	for assembly in assembly_list:
		# ------- Get HDF5 file handles list for each type of assembly and store in a dictionary.

		files  = {core.get_index(fname.filename.split('/')[-1]):core.get_files(core.get_directory(str(assembly)),fname.filename.split('/')[-1]) for fname in core.get_files(core.get_directory(str(assembly)))}
		index_redshift_dict  	= core.get_index_redshift_dict(assembly)
		kappa_co_assembly  		= dict()
		for index in files:
			coords  				= files[index]['Coordinates']*u.Mpc
			masses 					= files[index]['Mass']*u.M_sun
			velocities  			= files[index]['Velocity']*(u.km/u.s)
			net_velocities  		= np.linalg.norm(velocities,axis=1)
			rxv						= np.cross(coords,velocities)
			angular_momentum  		= get_angular_momentum(rxv,masses)
			total_angular_momentum 	= np.sum(angular_momentum,axis=0)
			zaxis  					= total_angular_momentum/np.linalg.norm(total_angular_momentum)
			angular_momentum_mask  	= np.dot(angular_momentum,zaxis)>0
			kinetic_energy  		= 0.5*np.sum(masses*net_velocities**2)
			rotational_kinetic_energy  	= 0.5*np.sum((masses*net_velocities**2)[angular_momentum_mask])
			kappa_co_assembly[index_redshift_dict[index]] 	= (rotational_kinetic_energy/kinetic_energy).value

		kappa_co[assembly] = kappa_co_assembly
		
		for index in files:
			files[index].close()

	for assembly in kappa_co:
		sns.lineplot(x=kappa_co[assembly].keys(),y=kappa_co[assembly].values())
		plt.show()

	print(kappa_co)
