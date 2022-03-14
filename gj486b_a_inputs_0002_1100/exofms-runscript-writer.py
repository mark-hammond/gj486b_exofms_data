#! /usr/bin/env python

'''
exofms-runscript-writer.py
MDH 25/06/19
Writes runscript, namelist, and diag_table and runs model
'''

import os
from write_constants import write_constants_fn
from exofmsRunUtils import write_runscript, write_namelist, write_diag_table, write_field_table
from run_params_class import run_params
import json


cwd = os.path.dirname(os.path.realpath(__file__))

for file in os.listdir(cwd):
    if file.startswith("inputs_"):
        inputs_f = open(file)

inputs = json.load(inputs_f)


# Search config file for output directory
with open('user_config', 'r') as f:
    for line in f:
        config = [x.strip() for x in line.split('=')]
        if config[0] == 'output_dir':            
            output_dir = config[1]
        if config[0] == 'fms_home':
            fms_home = config[1]
        if config[0] == 'restart_folder':
            restart_folder = config[1]

# Create model run object
exofms_run = run_params()

# --- Runscript parameters ---
# Run name
exofms_run.run_name = 'gj486b_a_'+inputs_f.name[:-5]
# Timestep in seconds
exofms_run.timestep_seconds = 30
# Run length in days
exofms_run.run_length_days = 1100
# Restart folder ("" for new run)
# exofms_run.restart_folder = restart_folder
exofms_run.restart_folder = ""
exofms_run.restart_days = 0


# --- Model parameters ---
# Cubed-sphere resolution
np = 24
# Number of vertical levels (see fv_eta for details)
exofms_run.fv_core_nml.npz = 32
# Number of x grid points
exofms_run.fv_core_nml.npx = np+1
# Number of y grid points
exofms_run.fv_core_nml.npy = np+1
# Number of CPUs
exofms_run.n_cpus = 6
# Cubed-sphere grid layout
exofms_run.fv_core_nml.layout = [1,1]
exofms_run.fv_core_nml.io_layout = [1,1]
# Print output summary frequency in hours
exofms_run.fv_core_nml.print_freq = 6
# Number of water-like substances (minimum 1)
exofms_run.fv_core_nml.nwat = 1

#  --- Constants ---
R = 8.31446261815324 
# Radius
radius = 8320526.0
# Angular velocity
omega = 0.00004956793
# Gravity
grav = 16.229474811
# Gas constant
rdgas = 8314.0 / float(inputs['mmw'])#287.0
# Lapse rate
kappa = 2./7.
# Mean molecular weight of air
wtmair = R*1000/rdgas

# -- Vertical levels ---
exofms_run.fv_eta_nml.uniform_z_spacing = False
exofms_run.fv_eta_nml.flexible_grid = False

# Atmosphere namelist
exofms_run.atmosphere_nml.do_dynamics = True
# Physical timestep split factor
exofms_run.atmosphere_nml.p_split = 1
# Dynamical timestep split factor (0 to let the model calculate this)
exofms_run.fv_core_nml.n_split = 0

#  --- Physics parameters ---
exofms_run.exo_phys_nml.relax_module = False
exofms_run.exo_phys_nml.do_vert_diff_bl = True
exofms_run.exo_phys_nml.do_simple_bl = False
exofms_run.exo_phys_nml.do_dry_adjustment = True
exofms_run.exo_phys_nml.do_dry_convection = False
exofms_run.exo_phys_nml.surface_on = True
exofms_run.exo_phys_nml.cp_surf = 1.e5
exofms_run.exo_phys_nml.tidally_locked = True


#  --- Initialisation parameters ---
# Type of initialisation (1 - isotherm, 2 - Guillot 2010, 3 - adiabat)
exofms_run.initialisation_nml.init_type = 3
# Surface pressure in Pa
exofms_run.initialisation_nml.surface_pressure = float(inputs['p0'])
# Surface temperature
exofms_run.initialisation_nml.ts = 600.0
# Stratosphere temperature
exofms_run.initialisation_nml.T_strat = 300.0


# --- Dry convection parameters ---
exofms_run.dry_convection_nml.tau = 100000.0
exofms_run.dry_convection_nml.gamma = 1.0

# -- Grey-gas radiation parameters ---
exofms_run.rad_coupler_nml.I0 = 55130.4
exofms_run.rad_coupler_nml.Tint = 0.0
exofms_run.rad_coupler_nml.rad_scheme = 'ts_twostream'
exofms_run.rad_coupler_nml.tau_IRe0 = float(inputs['taulw'])
exofms_run.rad_coupler_nml.tau_Ve0 = 0.0#float(inputs['tausw'])

# Order of damping
# Divergence damping (1: del-4, 2: del-6, 3: del-8)
# Vorticity diffusion (1: del-4, 2: del-6, 3: del-6)
exofms_run.fv_core_nml.nord = 1
# Divergence damping coefficient (~ 0.1 - 0.15)
exofms_run.fv_core_nml.d4_bg = 0.12
# Vorticity diffusion coefficient (~ 0.02 - 0.06)
exofms_run.fv_core_nml.vtdm4 = 0.03


# coefficient for external (barotropic) mode damping
exofms_run.fv_core_nml.d_ext = 0.02
# Number of layers from the top on which 2z energy-momentum-mass conserving filter is applied (recommed above 100mb)
exofms_run.fv_core_nml.n_sponge = 1
# Rayleigh drag profile center
exofms_run.fv_core_nml.rf_center = 3.e0
# Rayleigh drag ref pressure
exofms_run.fv_core_nml.p_ref = exofms_run.initialisation_nml.surface_pressure
# Rayleigh drag timescale
exofms_run.fv_core_nml.tau = 0.0


# -- Diffusive boundary layer parameters --
exofms_run.diffusivity_nml.depth_0 = 5000.0
exofms_run.diffusivity_nml.fixed_depth = True
exofms_run.diffusivity_nml.free_atm_diff = False
exofms_run.mixed_layer_nml.depth = 100.0

# -- Tracers --
# Leave sphum otherwise code crashes
exofms_run.field_table.tracers.append('sphum')
# Add any other required tracers below:
#exofms_run.field_table.tracers.append('cloud_water')

# -- Diagnostics ---
# Add diagnostic files
exofms_run.diag_table.lines.append('"atmos_static", -1,  "hours",  1, "days", "time",')
exofms_run.diag_table.lines.append('"atmos_average", -1,  "hours",  1, "days", "time",')
exofms_run.diag_table.lines.append('"atmos_daily",  2400  "hours",  1, "days", "time",')

# Add default diagnostic fields
exofms_run.diag_table.lines.append('"dynamics", "bk"         "bk"         "atmos_static", "all", .false.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "pk"         "pk"         "atmos_static", "all", .false.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "hyam"       "hyam"       "atmos_static", "all", .false.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "hybm"       "hybm"       "atmos_static", "all", .false.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "zsurf"      "zsurf"      "atmos_static", "all", .false.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "grid_lon",  "grid_lon",  "atmos_static", "all", .false.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "grid_lat",  "grid_lat",  "atmos_static", "all", .false.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "grid_lont", "grid_lont", "atmos_static", "all", .false.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "grid_latt", "grid_latt", "atmos_static", "all", .false.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "area",      "area",      "atmos_static", "all", .false.,  "none", 2,')

# Add diagnostic fields
exofms_run.diag_table.lines.append('"dynamics", "ps"         "ps"         "atmos_average", "all", .true.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "temp",      "temp",      "atmos_average",  "all", .true.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "ppt",      "ppt",      "atmos_average",  "all", .true.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "height",     "height",     "atmos_average",  "all", .true.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "ucomp",      "ucomp",      "atmos_average",  "all", .true.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "vcomp",      "vcomp",      "atmos_average",  "all", .true.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "omega",      "omega",      "atmos_average",  "all", .true.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "tsurf",      "tsurf",      "atmos_average",  "all", .true.,  "none", 2,')

exofms_run.diag_table.lines.append('"dynamics", "ps",     "ps",     "atmos_daily",  "all", .true.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "temp",     "temp",     "atmos_daily",  "all", .true.,  "none", 2,')
#exofms_run.diag_table.lines.append('"dynamics", "ppt",     "ppt",     "atmos_daily",  "all", .true.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "tsurf",     "tsurf",     "atmos_daily",  "all", .true.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "height",     "height",     "atmos_daily",  "all", .true.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "ucomp",     "ucomp",     "atmos_daily",  "all", .true.,  "none", 2,')
exofms_run.diag_table.lines.append('"dynamics", "vcomp",     "vcomp",     "atmos_daily",  "all", .true.,  "none", 2,')
#exofms_run.diag_table.lines.append('"dynamics", "omega",     "omega",     "atmos_daily",  "all", .true.,  "none", 2,')
exofms_run.diag_table.lines.append('"ts_short_char_bg", "olr",      "olr",      "atmos_daily",  "all", .true.,  "none", 2,')

# Each of the six faces has dimensions np*np
# Each face is then divided as per "layout"
# Note that the number of cells (layout[0]*layout[1]*6) 
# must be an integer multiple of n_cpus

# npx-1 must divide by layout[0]
# npy-1 must divide by layout[1]
exofms_run.resolutionXX = str(np)
exofms_run.resolutionCXX = 'C'+str(np)

if (exofms_run.fv_core_nml.npx-1)%exofms_run.fv_core_nml.layout[0]:
    print('Invalid grid')
    quit()
if (exofms_run.fv_core_nml.npy-1)%exofms_run.fv_core_nml.layout[1]:
    print('Invalid grid')
    quit()


write_constants_fn(radius,omega,grav,rdgas,kappa,wtmair)

# Write namelist
write_namelist(exofms_run)

# Write diag table
write_diag_table(exofms_run)

# Write field table
write_field_table(exofms_run)

# Write runscript
write_runscript(exofms_run)

# Run model
err = os.system('./run.sh ' + fms_home)

if (err !=0 ):
    exit('ERROR: Model compilation, run or horizontal regridding has failed')

# Print finished
print("Finished run")

# --- Move output data ---
experiment_output_dir = output_dir + '/' + exofms_run.run_name
os.system('mkdir -p ' + experiment_output_dir)

experiment_run_output_dir = output_dir + '/' + exofms_run.run_name + '/' + exofms_run.run_name + '_' + str(exofms_run.run_length_days+exofms_run.restart_days)

if os.path.isdir(experiment_run_output_dir) == True:
	os.system('rm -rf ' + experiment_run_output_dir)

os.system('mkdir ' + experiment_run_output_dir)
os.system('rm -rf ' + 'workdir/atmos*tile*.nc')
os.system('mv ' + 'workdir/atmos_* ' + experiment_run_output_dir + '/')
os.system('cp ' + 'exofms-runscript-writer.py ' + experiment_run_output_dir + '/')
os.system('mv ' + 'workdir/RESTART ' + experiment_run_output_dir + '/')

save_srcmods = True
if save_srcmods == True:
	os.system('cp -rf ' + 'srcmods ' + experiment_run_output_dir + '/')


print("Moved data")

# Clean run directory
#os.system('rm -rf workdir')

# Load modules necessary to run vertical regridding
module_home = os.getenv('MODULESHOME')
exec(open(module_home+'/init/python.py').read())
module('load', 'intel-compilers/2018upd4')
module('load', 'intel-mpi/2018upd4')
module('load', 'netcdf/netcdf-c-4.7.3')
module('load', 'netcdf/netcdf-fortran-4.5.2')

regrid_path = fms_home+'/exofms-utils/regrid_vertical'
print('regrid path: ', regrid_path )

# Alter any paths necessary in the vertical regridding
os.system('sed -i "s|^\s*interper.*|    interper = \'' + regrid_path +'/scripts/plevel.sh\'|g" '+regrid_path+'/scripts/plevel_fn.py')
os.system('sed -i "s|^executable=.*|executable=' + regrid_path +'/exec/plev.x|g" '+regrid_path+'/scripts/plevel.sh')

# Compile the vertical regridding
err = os.system(regrid_path+'/compile_plev_interpolation.sh ' + regrid_path)

if (err != 0):
    exit('ERROR: Compilation of vertical regridding failed')

# Run vertical regridding
err = os.system('python '+ fms_home + '/exofms-utils/regrid_vertical/scripts/interp.py ' + experiment_run_output_dir + '/ ' + str(exofms_run.radiation_nml.reference_slp))

if (err != 0):
    exit('ERROR: Vertical regridding failed')

print("Vertically interpolated")


inputs_f.close()

