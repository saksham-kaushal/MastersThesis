import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import seaborn as sns

import core

# =============================== User-defined functions ======================

def get_angular_momentum(rxv, mass):
	return np.transpose(np.multiply(mass, np.transpose(rxv)))

def thob_implementation(coords, velocities, masses):
	distances 		= np.linalg.norm(coords,axis=1)						# 114
	smomentums  	= np.cross(coords,velocities) 						# 126
	momentum  		= np.sum(masses[:,np.newaxis]*smomentums,axis=0) 	# 127
	Momentum 		= np.linalg.norm(momentum)							# 128
	zaxis 			= (momentum/Momentum)								# 130
	zheight 		= np.sum(zaxis*coords,axis=1)						# 131
	cylposition 	= coords-zheight[:,np.newaxis]*[zaxis]				# 132
	cyldistances 	= np.sqrt(distances**2-zheight**2)					# 133
	smomentumz 		= np.sum(zaxis*smomentums,axis=1)					# 134
	vrots 			= smomentumz/cyldistances							# 135
	Mvrot2 			= np.sum((masses*vrots**2)[vrots>0])				# 139
	kappa 			= Mvrot2 / np.sum(masses *					# 140
    	(np.linalg.norm(velocities,axis=1))**2)
	return kappa.value


# =============================== Main program ======================

if __name__ == '__main__':

	assembly_list 			= ['GM-Early','Organic','GM-Late']
	# assembly_list 			= ['Organic']
	distance_limit  		= 0.03*u.Mpc

	kappa_co 		= dict()
	kappa_co_thob  	= dict()

	for assembly in assembly_list:
		# ------- Get HDF5 file handles list for each type of assembly and store in a dictionary.

		files  = {core.get_index(
			fname.filename.split('/')[-1]):core.get_files(
			core.get_directory(
				str(assembly)),fname.filename.split('/')[-1]) 
			for fname in core.get_files(
				core.get_directory(
					str(assembly)))}

		index_redshift_dict  	= core.get_index_redshift_dict(assembly)
		kappa_co_assembly  		= dict()
		kappa_co_assembly_thob  = dict()
		for index in files:
			coords  				= files[index]['Coordinates']*u.Mpc
			masses 					= files[index]['Mass']*u.M_sun
			velocities  			= files[index]['Velocity']*(u.km/u.s)
			net_velocities  		= np.linalg.norm(velocities,axis=1)
			rxv						= np.cross(coords,velocities)
			angular_momenta  		= get_angular_momentum(rxv,masses)
			total_angular_momentum 	= np.sum(angular_momenta,axis=0)
			zaxis  					= total_angular_momentum/np.linalg.norm(total_angular_momentum)
			angular_momentum_mask  	= np.dot(angular_momenta,zaxis)
			kinetic_energy  		= 0.5*np.sum(masses*net_velocities**2)
			rotational_kinetic_energy  	= 0.5*np.sum((masses*net_velocities**2)[angular_momentum_mask>0])
			kappa_co_assembly[index_redshift_dict[index]] 	= (rotational_kinetic_energy/kinetic_energy).value
			kappa_co_assembly_thob[index_redshift_dict[index]] 	= thob_implementation(coords,velocities,masses)

		kappa_co[assembly] 		= kappa_co_assembly
		kappa_co_thob[assembly] = kappa_co_assembly_thob
		
		for index in files:
			files[index].close()

	for assembly in kappa_co:
		sns.lineplot(x=kappa_co[assembly].keys(),y=kappa_co[assembly].values())
		sns.lineplot(x=kappa_co_thob[assembly].keys(),y=kappa_co_thob[assembly].values())
		plt.show()

	print(kappa_co)
	print(kappa_co_thob)


