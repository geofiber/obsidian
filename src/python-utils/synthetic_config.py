#!/usr/bin/env python

"""
RS 2018/04/12:  Config Settings for Synthetic Dataset

This problem is hard and we need synthetic data to understand what all
is going on.  So here's a generator for an easy synthetic configuration.
"""

import copy
import numpy as np
import pandas as pd
from config_writer import write_config

# ============================================================================
#                USER-CONFIGURABLE OPTIONS WITHIN THIS BLOCK
# ============================================================================

# First we need to define a configuration, which is going to consist of
# latitude and longitude boundaries of the area, as well as a list of
# (name, min, max) for different layers.  We have to be careful about the
# order in which these layers are specified because the configuration
# fixes a strict order in which those layers can appear in the model,
# even if some layers are absent at certain locations.  Example below.

_datadir = "../../../data/"
source_petrofn = _datadir + "Petrophysics_Gascoyne.csv"
config_petrofn = _datadir + "Petrophysics_Synthetic.csv"
gravdata_fn = _datadir + "gravity_400m_Synthetic.txt"
magdata_fn = _datadir + "mag_TMI_Synthetic.txt"
fieldsrc_fn = _datadir + "Formation_data_Gascoyne.csv"
fielddata_fn = _datadir + "Formation_data_Synthetic.csv"

# Column format for sensor files
sensor_colnames = ['id', 'val', 'lat', 'lng']

# Ambient magnetic field from the IGRF model, which can be queried at
#     https://www.ngdc.noaa.gov/geomag-web/#igrfwmm
H_IGRF = np.array([1.419e+2, 2.8739e+4, -4.62667e+4])

config_layers = pd.DataFrame(
        # layer name, layer type, (min, max) depth in m, (nx, ny) control pts
        [('Layer A', 'normal',  -1.0e+0, 1.0e+0, 1, 1),
         ('Layer B',  'normal', -3.0e+3, 3.0e+3, 1, 2)],
        columns=['name','type','zmin','zmax','nx','ny'])

config_params = { 'lng': 116.10, 'lat': -24.85, 'L': 2.0e+4,
                  'maxdepth': 0.1e+4, 'layers': config_layers,
                  'H_IGRF': H_IGRF, }

# Rock properties
layer_rockprops = pd.DataFrame(
        # layer name, mu_rho, mu_chi_M
        [('Layer A', 'Durlacher Supersuite', 2.6, 0.5),
         ('Layer B', 'Halfway Gneiss',       3.0, 2.5)],
        columns=['name', 'subname', 'rho_mean', 'chiM_mean'])

# ============================================================================
# REST OF CODE:  calculates rock priors, reads in data, & writes configuration
# ============================================================================

def generate_rock_data(layer_rockprops, source_petrofn, config_petrofn):
    """
    Generate some synthetic rock data based on well-defined properties.
    Quick way to do this:  rip off existing petrophysics data, keep rock
    metadata, and overwrite existing rock properties w/synthetic data.
    """
    layers_to_cat = [ ]
    rockdata = pd.read_csv(source_petrofn)
    psig, Xsig = 0.05, 0.25
    for i in range(len(layer_rockprops)):
        name, subname, pmu, Xmu = layer_rockprops.iloc[i,:]
        layerdata = pd.DataFrame(rockdata[rockdata.FormationName==subname])
        layerdata.loc[:,'FormationName'] = name
        layerdata.loc[:,'UnitName'] = name
        layerdata.loc[:,'Magnetic_susceptibility'] = 10**np.random.normal(
                Xmu, Xsig, size=len(layerdata))
        layerdata.loc[:,'Density_g_cm-3'] = np.random.normal(
                pmu, psig, size=len(layerdata))
        layers_to_cat.append(layerdata)
    rockdata_synth = pd.concat(layers_to_cat, axis=0)
    with open(config_petrofn, 'w') as outfile:
        outfile.write(rockdata_synth.to_csv(index=False))

def generate_field_data(layer_rockprops, fielddata_fn, fieldsrc_fn):
    """
    Generate some synthetic field data based on well-defined properties.
    Quick way to do this:  rip off existing petrophysics data, keep rock
    metadata, and overwrite existing rock properties w/synthetic data.
    """
    layers_to_cat = [ ]
    fieldsrc = convert_ground_truth(fieldsrc_fn)
    for i in range(len(layer_rockprops)):
        name, subname, pmu, Xmu = layer_rockprops.iloc[i,:]
        layerdata = pd.DataFrame(fieldsrc[fieldsrc.form_name==subname])
        layerdata.loc[:,'form_name'] = name
        layerdata.loc[:,'unit_name'] = name
        layers_to_cat.append(layerdata)
    fielddata_synth = pd.concat(layers_to_cat, axis=0)
    with open(fielddata_fn, 'w') as outfile:
        outfile.write(fielddata_synth.to_csv(index=False))
    return fielddata_synth

def compute_rock_priors(rockfname):
    """
    Compute distributions of rock properties for geological layers from
    petrophysical data provided by Hugo.  The prior has 9 components:
        Density in kg/m^3
        log(Magnetic susceptibility)
        Thermal conductivity
        Thermal productivity
        log(Resistivity x-component)
        log(Resistivity y-component)
        log(Resistivity z-component)
        P-wave velocity

    :param rockfname:  name of CSV file with rock measurements from Hugo
    :returns:  two dicts, rockprior_mu and rockprior_cov, containing the
        sample mean and covariance of the rock property measurements for
        each rock layer (indexed by the rock layer's name)
    """

    # Read in the petrophysics rock data Hugo gave us.
    rockdata = pd.read_csv(rockfname)
    rockpriormu = { }
    rockpriorcov = { }
    # Default janky prior just in case we don't have enough physical samples.
    # The first two terms are derived from the petrophysics data aggregated
    # over all layers.  The other terms are made up because they have to be
    # filled, but if we want to use other sensors we need better priors.
    rockprops = ['Density_kgm3', 'log_MagSusc',
                 'ThermConduct', 'ThermProduct',
                 'log_Resist_x', 'log_Resist_y', 'log_Resist_z',
                 'Resist_Phase', 'PWaveVelocity']
    rockpriormu_def = pd.Series(
            [ 2.7e+3, 1.5, 2.0, 2e-6, 0.0, 0.0, 0.0, 0.0, 4000 ],
            index=rockprops)
    rockpriorcov_def = pd.DataFrame(np.diag(
            [ 5.0e+4, 0.5, 0.25, 1e-12, 1.0, 1.0, 1.0, 1.0, 2.5e+5 ]),
            index=rockprops, columns=rockprops)
    
    # Go through all the formations we have data for.
    forms = np.unique(rockdata.FormationName)
    for i, f in enumerate(forms):
        df = rockdata[rockdata.FormationName==f]
        rockpriormu[f] = copy.deepcopy(rockpriormu_def)
        rockpriorcov[f] = copy.deepcopy(rockpriorcov_def)
        # If we have enough data points (I say 5, somewhat arbitrarily),
        # construct an empirical Gaussian prior from Hugo's data.
        if len(df) >= 5:
            # Grab last two columns -- susceptibility and density
            dc = np.zeros(shape=(df.shape[0], 2))
            # Convert magnetic susceptibility to log scale
            dc[:,0] = 1000.0*df.loc[:,'Density_g_cm-3']
            # Convert density to kg/m^3
            dc[:,1] = np.log10(df.loc[:,'Magnetic_susceptibility'])
            # Store sample covariance in hash
            rockpriormu[f].iloc[:2] = np.mean(dc, axis=0)
            rockpriorcov[f].iloc[:2,:2] = np.cov(dc.T)
    return rockpriormu, rockpriorcov

def convert_ground_truth(fname):
    """
    Reads Hugo's formation data file and converts to a standard format.
    :param fname:  file with formation boundary information
    """
    fieldcols = ['id', 'site_id', 'lat', 'lng', 'form_name',
                 'unit_name', 'geochron', 'sample_id', 'age', 'age_err']
    fieldtypes = [str, str, float, float, str, str, str, str, float, float]
    dtype = { fc: ft for fc, ft in zip(fieldcols, fieldtypes) }
    fielddata = pd.read_csv(fname, names=fieldcols,
                            dtype=dtype, na_values=['-'], skiprows=1)

    # Calculate a scalar index to represent the layer boundaries, based on
    # what we've got in our data files so far.  Add as 'val' column to data.
    fieldval_lookup = { fn: i for i, fn in enumerate(config_layers.name) }
    fieldval = pd.Series([fieldval_lookup.get(fn, -1) for fn in fielddata.form_name])
    fielddata = pd.concat([fielddata, fieldval], axis=1)
    fielddata.rename(columns={0:'val'}, inplace=True)

    return fielddata

def main():
    """
    The main config-writer routine
    """

    # Generate rock priors
    generate_rock_data(layer_rockprops, source_petrofn, config_petrofn)
    rockpriormu, rockpriorcov = compute_rock_priors(config_petrofn)

    # Read prospector NPZ
    prospect = np.load('synthetic.npz')
    gravsynth = prospect['gravReadings'][0]
    gravsynth -= np.mean(gravsynth)
    # gravsynth += 0.05*np.std(gravsynth)*np.random.normal(size=gravsynth.shape)
    with open('gravsynth.csv', 'w') as csvfile:
        csvfile.write(pd.DataFrame(gravsynth).to_csv(index=False, header=False))
    magsynth = prospect['magReadings'][0]
    magsynth -= np.mean(magsynth)
    # magsynth += 0.05*np.std(magsynth)*np.random.normal(size=magsynth.shape)
    with open('magsynth.csv', 'w') as csvfile:
        csvfile.write(pd.DataFrame(magsynth).to_csv(index=False, header=False))

    # Read in grav data and convert to the standard format.
    # For Hugo's data, divide by 10 to convert from um/s^2 to mgal.
    gravdata = pd.read_csv(gravdata_fn, names=sensor_colnames,
                           dtype=float, skiprows=1)
    gravdata.val /= 10.0

    # Read in mag data and convert to the standard format.
    magdata = pd.read_csv(magdata_fn, names=sensor_colnames,
                          dtype=float, skiprows=1)

    # Read in field observation data and convert to the standard format.
    fielddata = generate_field_data(layer_rockprops, fielddata_fn, fieldsrc_fn)

    return

    # Fill remaining fields in config_params
    config_params.update({ 'grav_data': gravdata,
                           'mag_data': magdata,
                           'field_data': fielddata,
                           'rockpriormu': rockpriormu,
                           'rockpriorcov': rockpriorcov, })

    write_config(**config_params)

if __name__ == "__main__":
    main()
