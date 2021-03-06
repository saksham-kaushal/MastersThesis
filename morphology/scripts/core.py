import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import h5py
import os
import re
import inspect
from astropy.cosmology import FlatLambdaCDM, z_at_value


def get_directory(directory):
	'''
	Sets up directory structure. Returns directory path for root, data, scripts, plots, etc.
	Parameters :
	directory - Identifier for the directory label, eg. root, data, scripts, plots, GM-Early_data, Organic_data, GM-Late_data.
	'''
	dir_dict 						= dict()
	dir_dict['scripts_dir'] 		= os.path.dirname(os.path.abspath(__file__))
	dir_dict['root_dir'] 			= os.path.dirname(dir_dict['scripts_dir'])
	dir_dict['data_dir'] 			= os.path.join(dir_dict['root_dir'], 'data')
	dir_dict['GM-Early_dir']		= os.path.join(dir_dict['data_dir'], 'GM-Early')
	dir_dict['Organic_dir']			= os.path.join(dir_dict['data_dir'], 'Organic')
	dir_dict['GM-Late_dir']			= os.path.join(dir_dict['data_dir'], 'GM-Late')
	dir_dict['plots_dir']			= os.path.join(dir_dict['root_dir'], 'plots')
	dir_dict['zavala2l_dir']		= os.path.join(dir_dict['data_dir'], 'zavala_fig2l')
	dir_dict['baryon_cdm_data_dir']	= os.path.join(dir_dict['data_dir'], 'baryon_cdm_data')
	dir_dict['GM-Early_baryon_cdm_data_dir']	= os.path.join(dir_dict['baryon_cdm_data_dir'], 'GM-Early')
	dir_dict['Organic_baryon_cdm_data_dir']		= os.path.join(dir_dict['baryon_cdm_data_dir'], 'Organic')
	dir_dict['GM-Late_baryon_cdm_data_dir']		= os.path.join(dir_dict['baryon_cdm_data_dir'], 'GM-Late')
	
	if directory == 'plotter':							# If plotter function calls this function, return the plotting subdirectory in the "plots" directory.
		frame 				= inspect.stack()[-1]
		plot_subdir_path	= frame[0].f_code.co_filename
		plot_subdir 		= os.path.splitext(os.path.basename(plot_subdir_path))[0]
		return os.path.join(dir_dict['plots_dir'], str(plot_subdir))
	
	return dir_dict[directory+'_dir']


def get_files(directory,fname=None,mode='r'):
	'''
	Returns handle(s) for hdf5 file(s). If no filename is mentioned, handles for all files in the directory are returned; for a filename or a list of filenames corresponding file handle or a list of file handles is returned.
	Parameters	:
	directory 	- path to the directory containing hdf5 files.
	fname 		- name of a single or a list of hdf5 files to read. Defaults to None, so that all files in the directory are read. 
	mode 		- mode for opening file. Defaults to read mode.
	'''
	# assert isinstance(fname,list), 'Pass file names as list.'
	if fname is None:
		return [h5py.File(os.path.join(directory,fname),mode) for fname in os.listdir(directory)]
	elif isinstance(fname,list):
		return [h5py.File(os.path.join(directory,fname_handle),mode) for fname_handle in fname]
	elif isinstance(fname,str):
		return h5py.File(os.path.join(directory,fname),mode)

		
def get_redshift(fname):
	'''
	Returns redshift value for an hdf5 file. Uses regex to find redshift value from filename.
	Parameters	:
	fname 	- Open handle for an hdf5 file.
	'''

	return float(re.compile(r'z\d+p\d+').findall(str(fname))[0].replace('z','').replace('p','.'))


def get_index(fname):
	'''
	Returns index value for an hdf5 file. Uses regex to find redshift value from filename.
	Parameters	:
	fname 	- Open handle for an hdf5 file.
	'''
	return int(re.compile(r'_\d+').findall(fname)[0].replace('_',''))


def get_index_redshift_dict(assembly='GM-Early'):
	'''
	Returns a dictionary of index and redshift values from filenames for a given assembly mode.
	Parameters	:
	assembly 	- Assembly mode. Defaults to the GM-Early mode.
	'''
	assembly_dir  	= get_directory(assembly)
	index_calibration_file_list 	= [fname for fname in os.listdir(assembly_dir)]
	return dict(zip(map(get_index,index_calibration_file_list),map(get_redshift,index_calibration_file_list)))


def get_time_from_redshift(redshift,cosmology=None):
	'''
	Returns cosmological time since big-bang from redshift value for a flat lambda-CDM cosmology utilising parameters used by EAGLE simulation suite, unless specified.
	Parameters	:
	assembly 	- Assembly mode. Defaults to the GM-Early mode.
	'''
	if cosmology is None:
		cosmology 	= FlatLambdaCDM(100.*0.6777,Om0=0.307,Ob0=0.04825)
	return cosmology.age(redshift) 

def get_plot_axes_titles():
	'''
	Returns a dict containing plot axis titles with standard units for different physical parameter.
	'''
	title  									= dict()
	title['redshift']						= 'Redshift'
	title['net_specific_angular_momentum']	= 'Specific angular momentum (Mpc km s$^{-1}$)'
	title['expansion_factor']				= 'Expansion factor'
	title['net_angular_momentum']			= 'Angular momentum (Mpc M$_{\odot}$ km s$^{-1}$)'
	title['total_kinetic_energy']			= 'Kinetic energy (M$_{\odot}$ km$^{2}$ s$^{-2}$)'
	title['time'] 							= 'Time (Gyr)'
	title['masses']							= 'Mass (M$_{\odot}$)'
	title['median_radius']					= 'Median radius (Mpc)'
	return title
	
def prepare_plot(context='paper',theme='dark',font_scale=1,colours=None,rc_kwparams=dict()):
	'''
	Set seaborn styling for plots. Use print(sns.axes_style()) for getting a complete list of style attributes
	Parameters	:	
	context		- seaborn plot context. Defaults to paper.
	theme 		- seaborn plot theme. Defaults to dark theme.
	font_scale 	- font scaling for plot text. Defaults to 1.
	colours 	- list of plot marker colours. Defaults to None changed to magenta-green-blue later in the function.
	rc_kwparams	- dictionary of keyword parameters to be passed as matplotlib rc paramerters (which are essentially the runtime configuration settings).
	'''
	rc_params	= {
		'xtick.bottom': True, 
		'ytick.left': True, 
		'axes.spines.left': True, 
		'axes.spines.bottom': True
		}
	rc_params.update(rc_kwparams)			# Update parameters by adding parameters passed to function
	sns.set(style=theme,font='TeX Gyre Pagella Math')
	sns.set_context(context,font_scale=font_scale)
	sns.set_style(rc_params)
	if colours is None:
		colours  	= [
						'#772e51', 		# Red/Magenta (lighter alternative '#923563')
						'#44b677',		# Green
						'#66b6d2'		# Blue (alternatives '#88ccee',#94c8e0', '#51abcb')
					]
	sns.set_palette(sns.color_palette(colours))
	return None
	
def plot_or_not(show,plot_name=None,suffix='',dpi=480,ftype='png',bbox_inches='tight'):
	'''
	Function to switch between show and save methods for plots.
	Parameters	:	
	show		- Shows an 'interactive' plot if True, else saves the plot if False. If set to None, does nothing.
	plot_name	- If set, plot is saved using plot_name if show is set to False. If show is set to False, yet plot_name is not provided, a random number is generated for filename.
	dpi 		- dots per inch in the saved plot, default 480.
	ftype 		- file type of image to save, default .png.
	bbox_inches	- bounding box layout when saving a figure, default tight layout.
	'''
	if show == True :
		plt.show()
	elif show == False :
		if plot_name == None :
			plot_name = np.random.randint(10000,99999)
		plt.savefig(os.path.join(get_directory('plotter'),plot_name+str(suffix)+'.'+str(ftype)),
			dpi=dpi,bbox_inches=bbox_inches)
	elif show == None :
		pass 
