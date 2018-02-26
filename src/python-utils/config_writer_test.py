#!/usr/bin/env python

"""
RS 2018/01/24:  Script to write configuration files for Obsidian

This is a test suite for the configuration writing script.  That code is
long and elaborate enough that I feel like it needs some tests, even if
I don't think we're going to be messing with it a whole lot.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from config_writer import config_params, config_layers, config_fname
from config_writer import compute_rock_priors, form_csv_fnames, write_config
from config_writer import write_sensor_data, display_ground_truth

def test_compute_rock_priors():
    """
    Unit test of compute_rock_priors.
    """
    mu, cov = compute_rock_priors(config_petrofn)
    # If we got this far we must have finished, so assume we at least have
    # a hash full of means and covariances; let's see what it looks like
    for f in sorted(mu.keys()):
        print "Priors for {}:".format(f)
        print mu[f]
        print cov[f]

def test_form_csv_names():
    """
    Unit test of form_csv_names.
    """
    for prop in ['CtrlMean', 'CtrlMins', 'CtrlMaxs', 'CtrlMask',
                 'RockMean', 'RockMins', 'RockMaxs', 'RockMask']:
         print form_csv_fnames(config_layers, prop)

def test_write_sensor_data():
    """
    Unit test of write_sensor_data.
    """
    gravdata = pd.read_csv(
            "../../data/gravity_400m_Gascoyne.txt",
            names=['id','val','lat','lng'], dtype=float, skiprows=1)
    kwargs = { 'lng': config_params['lng'],
               'lat': config_params['lat'],
               'L':   config_params['L'], }
    write_sensor_data(gravdata, stag='grav', **kwargs)

def test_write_config():
    """
    Unit test of write_config.
    """
    write_config(**config_params)

def test_display_ground_truth():
    """
    Unit test of display_ground_truth.
    """
    display_ground_truth(
            config_params['lng'], config_params['lat'], config_params['L'])


if __name__ == "__main__":
    # test_compute_rock_priors()
    # test_form_csv_names()
    # test_write_sensor_data()
    test_write_config()
    # test_display_ground_truth()
