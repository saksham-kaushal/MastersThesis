import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import seaborn as sns

from astropy.cosmology import FlatLambdaCDM, z_at_value
cosmology 	= FlatLambdaCDM(100.*0.6777,Om0=0.307,Ob0=0.04825)

import core
import kappa_co as kco

# =============================== User-defined functions ======================

def get_specific_angular_momentum(coords, velocities, masses):
	'''
	Returns total angular momentum per unit mass.
	'''
	rxv							= np.cross(coords,velocities)
	angular_momenta  			= kco.get_angular_momentum(rxv,masses)
	total_angular_momentum 		= np.sum(angular_momenta,axis=0)
	specific_angular_momentum 	= total_angular_momentum/np.sum(masses)
	return np.linalg.norm(specific_angular_momentum)

def plot_kappa_specific_l(kappa_co,specific_angular_momentum,show=True,plot_name='kappa_specific_L'):
	'''
	Plots kappa versus speific angular momentum.
	'''
	core.prepare_plot(context='talk',font_scale=.9)
	fig  	= plt.figure()
	axis 	= fig.add_subplot(1,1,1)
	for assembly in kappa_co:
		assembly_plot_dict = dict()
		for key in kappa_co[assembly]:
			assembly_plot_dict[kappa_co[assembly][key]] = specific_angular_momentum[assembly][key].value
		print(assembly_plot_dict)
		x,y = zip(*sorted(assembly_plot_dict.items()))
		axis.plot(x,y,label=assembly)
	axis.legend(title='Assembly')
	axis.set_xlabel('$\kappa_{CO}$')
	axis.set_ylabel(core.get_plot_axes_titles()['net_specific_angular_momentum'])
	core.plot_or_not(show=show,plot_name=plot_name)

# =============================== Main program ======================

if __name__ == '__main__':

	assembly_list 			= ['GM-Early','Organic','GM-Late']
	# assembly_list 			= ['Organic']

	kappa_co_thob_et_al 		= dict()
	specific_angular_momentum 	= dict()

	for assembly in assembly_list:
		# ------- Get HDF5 file handles list for each type of assembly and store in a dictionary.

		files  = {core.get_index(
			fname.filename.split('/')[-1]):core.get_files(
			core.get_directory(
				str(assembly)),fname.filename.split('/')[-1]) 
			for fname in core.get_files(
				core.get_directory(
					str(assembly)))}

		index_redshift_dict  				= core.get_index_redshift_dict(assembly)
		kappa_co_assembly_thob_et_al 		= dict()
		specific_angular_momentum_assembly 	= dict()
		for index in files:
			coords  					= files[index]['Coordinates']*u.Mpc
			masses 						= files[index]['Mass']*u.M_sun
			velocities  				= files[index]['Velocity']*(u.km/u.s)
			kappa_co_assembly_thob_et_al[core.get_time_from_redshift(index_redshift_dict[index]).value] = kco.get_kappa_thob_et_al(coords,velocities, masses)
			specific_angular_momentum_assembly[core.get_time_from_redshift(index_redshift_dict[index]).value] = get_specific_angular_momentum(coords,velocities,masses)

		kappa_co_thob_et_al[assembly] 		= kappa_co_assembly_thob_et_al
		specific_angular_momentum[assembly] = specific_angular_momentum_assembly

		for index in files:
			files[index].close()

	plot_kappa_specific_l(kappa_co_thob_et_al,specific_angular_momentum,show=False)