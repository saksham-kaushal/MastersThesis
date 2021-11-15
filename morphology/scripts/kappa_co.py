import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import seaborn as sns

from astropy.cosmology import FlatLambdaCDM, z_at_value
cosmology 	= FlatLambdaCDM(100.*0.6777,Om0=0.307,Ob0=0.04825)

import core

# =============================== User-defined functions ======================

def get_eagle_cosmology():
	'''
	Returns flat lambda-CDM cosmology for astropy, using parameters used by EAGLE simulation suite. 
	'''
	return FlatLambdaCDM(100.*0.6777,Om0=0.307,Ob0=0.04825)

def get_angular_momentum(rxv, mass):
	'''
	Returns angular momentum three-vectors from mass and (cross-product of) radius and velocity.
	'''
	return np.transpose(np.multiply(mass, np.transpose(rxv)))

def get_orthonormal_basis(angular_momentum_vector):
	'''
	Returns a set of orthonormal basis vectors generated using a reference initial vector.
	'''
	k = angular_momentum_vector

	try:
		nonzero_k = np.nonzero(k)[0][0]
	except IndexError:
		raise IndexError('Received a null angular momentum vector')
	if nonzero_k == 0:
		i = np.array([-k[1]/k[0],1,0])
	else:
		i = np.array([1,0,0])
	
	j = np.cross(k,i)
	basis_vectors = np.vstack((i,j,k))
	return basis_vectors/np.linalg.norm(basis_vectors,axis=1)[:,np.newaxis]

def coordinates_transform(galaxy_coords,unit_3dvectors):
	'''
	Transforms three vectors to a new orthonormal frame of reference.
	'''
	transforms = list()
	for coord in galaxy_coords:
		transforms.append(np.dot(unit_3dvectors,coord))
	return np.array(transforms)

def get_kappa(coords, velocities, masses):
	'''
	Untested implementation to obtain value of kappa computed using two-dimensional plane projection method and coordinate transform of velocity values.
	'''
	net_velocities  		= np.linalg.norm(velocities,axis=1)
	rxv						= np.cross(coords,velocities)
	angular_momenta  		= get_angular_momentum(rxv,masses)
	total_angular_momentum 	= np.sum(angular_momenta,axis=0)
	basis_vectors 			= get_orthonormal_basis(total_angular_momentum.value)
	velocities_transform  	= coordinates_transform(velocities.value,basis_vectors)*(velocities.unit)
	velocities_transform[:,2]	= 0
	net_velocities_transform	= np.linalg.norm(velocities_transform,axis=1)
	angular_momentum_mask  	= np.dot(angular_momenta,basis_vectors[2])
	kinetic_energy  		= np.sum(masses*net_velocities**2)
	rotational_kinetic_energy  	= np.sum((masses*net_velocities_transform**2)[angular_momentum_mask>0])
	kappa 					= rotational_kinetic_energy/kinetic_energy
	return kappa.value

def get_kappa_thob_et_al(coords, velocities, masses):
	'''
	Returns the value of kappa computed using two-dimensional plane projection method employing direct use of the kappa equation by Correa et al and Thob et al.
	'''
	net_velocities  		= np.linalg.norm(velocities,axis=1)
	rxv						= np.cross(coords,velocities)
	angular_momenta  		= get_angular_momentum(rxv,masses)
	total_angular_momentum 	= np.sum(angular_momenta,axis=0)
	basis_vectors 			= get_orthonormal_basis(total_angular_momentum.value)
	coords_transform  		= coordinates_transform(coords.value,basis_vectors)*(coords.unit)
	coords_transform[:,2] 	= 0
	radius 					= np.linalg.norm(coords_transform,axis=1)
	angular_momentum_mask  	= np.dot(angular_momenta,basis_vectors[2])
	kinetic_energy  		= np.sum(masses*net_velocities**2)
	rotational_kinetic_energy  	= np.sum((np.dot(angular_momenta,basis_vectors[2])**2/(masses*radius**2))[angular_momentum_mask>0])
	kappa 					= rotational_kinetic_energy/kinetic_energy
	return kappa.value

def get_kappa_thob_implementation(coords, velocities, masses):
	'''
	Returns the value of kappa computed using the publicly available codes by Thob et al.
	'''
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
	kappa 			= (Mvrot2/2) / np.sum((masses *
    				  (np.linalg.norm(velocities,axis=1))**2)/2)		# 140
	return kappa.value


def kappa_plotter(kappa_co, show=True,plot_name='kappa'):
	'''
	Plots values of kappa vs time since big-bang. A secondary x-axis of redshift values is also used.
	'''
	core.prepare_plot(context='talk',font_scale=.9)
	fig  	= plt.figure()
	axis 	= fig.add_subplot(1,1,1)
	secondary_axis 		= axis.twiny()
	redshift_markers 	= [7.0,3.0,2.0,1.0,0.5,0.2,0.0]
	corresponding_times	= core.get_time_from_redshift(redshift_markers).value
	secondary_axis.set_xlim(axis.get_xlim())
	secondary_axis.set_xticks(corresponding_times)
	secondary_axis.set_xticklabels('{:02}'.format(redshift) for redshift in redshift_markers)
	secondary_axis.set_xlabel('Redshift')
	axis.set_xlabel('Time (Gyr)')
	axis.set_ylabel('$\kappa_{CO}^{}$')
	for assembly in kappa_co:
		sns.lineplot(x=kappa_co[assembly].keys(),y=kappa_co[assembly].values(),ax=axis)
	axis.legend(labels=kappa_co.keys())
	axis.plot(np.linspace(axis.get_xlim()[0],axis.get_xlim()[1],20),[0.4]*20,'--',color='black')
	core.plot_or_not(show=show, plot_name=plot_name)


# =============================== Main program ======================

if __name__ == '__main__':

	assembly_list 			= ['GM-Early','Organic','GM-Late']
	# assembly_list 			= ['Organic']

	kappa_co 				= dict()
	kappa_co_thob  			= dict()
	kappa_co_thob_et_al 	= dict()

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
		kappa_co_assembly_thob_et_al 	= dict()
		for index in files:
			coords  				= files[index]['Coordinates']*u.Mpc
			masses 					= files[index]['Mass']*u.M_sun
			velocities  			= files[index]['Velocity']*(u.km/u.s)
			# kappa_co_assembly[core.get_time_from_redshift(index_redshift_dict[index]).value] = get_kappa(coords,velocities, masses)
			kappa_co_assembly_thob_et_al[core.get_time_from_redshift(index_redshift_dict[index]).value] = get_kappa_thob_et_al(coords,velocities, masses)
			kappa_co_assembly_thob[core.get_time_from_redshift(index_redshift_dict[index]).value] = get_kappa_thob_implementation(coords,velocities,masses)

		# kappa_co[assembly] 				= kappa_co_assembly
		kappa_co_thob[assembly] 		= kappa_co_assembly_thob
		kappa_co_thob_et_al[assembly] 	= kappa_co_assembly_thob_et_al


		for index in files:
			files[index].close()
	
	# print(kappa_co_thob)
	# print(kappa_co)

	kappa_plotter(kappa_co_thob,show=False,plot_name='kappa_co_thob')
	kappa_plotter(kappa_co_thob_et_al,show=False,plot_name='kappa_co_thob_et_al')
	# kappa_plotter(kappa_co,show=False,plot_name='kappa_co')