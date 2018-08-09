#!/usr/bin/env python

"""
RS 2018/02/23:  Digesting Obsidian Output

Plots we want to see:
    Trace plots of key parameters like layer depths, rock properties
    Histograms of residuals between forward model values and data
    Contour maps of residuals between forward models and data
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import GPy

pickaxe_npz = ('/Users/rscalzo/Desktop/Projects/Professional/Geoscience'
               '/Formation Boundaries/sandbox/moomba-dk-run06.npz')
pickaxe_npz = ('/Users/rscalzo/Desktop/Projects/Professional/Geoscience'
               '/Formation Boundaries/sandbox/gascoyne_v4-thin1000.npz')
mason_npz   = ('/Users/rscalzo/Desktop/Projects/Professional/Geoscience'
               '/Formation Boundaries/sandbox/moomba-dk-run01.npz')


def plot_sensor(sensors, readings, chain, sample=None,
                units='unknown units', show=True, energy=None):
    """
    Plots the value of a sensor against the forward model:
        Contour plots of real data and forward models
        Contour plots of mean model residuals across the chain
        Histograms of model residuals across the chain
    :param sensors:  np.array of shape (N, 3) with the physical
        coordinates (x,y,z) of the sensors in meters; assumed to have a
        regular raster-like structure where N factors into (Nx, Ny)
    :param readings:  np.array of shape (N, ) with the observed
        sensor readings for the real dataset
    :param chain:  np.array of shape (M, N) with the synthetic
        sensor readings for each forward model
    :param sample:  int index of sample to grab from chain
        (defaults to None, which averages over the chain)
    :param units:  str describing the units of sensor readings
    :param show:  call plt.show()?  default True; set to False if you're
        making a multi-panel plot or want to save fig in calling routine
    :param energy:  optional log posterior estimate(s)
    """
    x, y, z = sensors.T
    d = readings - readings.mean()
    if energy is not None:    # HACK:  take MAP estimate
        sample = np.argmin(energy)
        print "Taking MAP estimate, energy =", energy[sample]
    if sample is None:
        print "Averaged fwd models over chain of shape", chain.shape
        f = chain.mean(axis=0) - chain.mean()
    elif not isinstance(sample, int):
        print "ERROR:  sample = {} is of type {}, not int".format(
            sample, sample.__class__)
        return
    elif abs(sample) > len(chain):
        print "ERROR:  sample = {} not in range ({}, {})".format(
                sample, -len(chain), len(chain))
        return
    else:
        print "Picking sample", sample, "from chain of shape", chain.shape
        f = chain[sample] - chain[sample].mean()

    # Contour map of residuals in f
    fig = plt.figure(figsize=(6,7))
    ax1 = plt.subplot2grid((3,1), (0,0), rowspan=1)
    plt.tricontourf(x, y, d, alpha=0.5, label='Data')
    plt.colorbar()
    plt.tricontour(x, y, f, colors='k', label='Fwd Model')
    plt.xlabel("Eastings (m)")
    plt.ylabel("Northings (m)")
    plt.legend(loc='upper right')

    # Contour map of residuals in f
    ax1 = plt.subplot2grid((3,1), (1,0), rowspan=1)
    plt.tricontourf(x, y, d-f, alpha=0.95, label='Data', cmap='coolwarm')
    plt.colorbar()
    plt.tricontour(x, y, f, colors='k', label='Fwd Model')
    plt.xlabel("Eastings (m)")
    plt.ylabel("Northings (m)")
    plt.legend(loc='upper right')

    # Histogram of residuals in f
    ax2 = plt.subplot2grid((3,1), (2,0), rowspan=1)
    plt.hist(d-f, bins=20)
    plt.xlabel("Data $-$ Model ({})".format(units))

    # Show
    plt.subplots_adjust(left=0.12, bottom=0.08,
                        right=0.90, top=0.92,
                        wspace=0.20, hspace=0.40)
    if show:
        plt.show()

def display_ground_truth(labels, spherical=True, show=True):
    """
    Displays geological ground-truth labels in a given area.  Put here
    until I find a better home for it.
    :param x:  x-coordinate of centre of extracted area
    :param y:  y-coordinate of centre of extracted area
    :param L:  length of side of (square) modeled area in metres
    :param labels:  pd.DataFrame with columns ['lat', 'lng', 'val']
    :param spherical:  bool, True if (x, y) = (lat, lng) in degrees
    :param show:  call plt.show()?  default True; set to False if you're
        making a multi-panel plot or want to save fig in calling routine
    """
    # Plot the ground truth
    x, y, v = labels.x, labels.y, labels.val
    for f in np.unique(v.values):
        idx = (v == f)
        plt.plot(x[idx], y[idx], ls='None', marker='o', ms=3, label=f)
    plt.legend(loc='upper left')
    # plt.title("${:.1f} \\times {:.1f}$ km$^2$ area centered on "
    #           "lng = ${:.3f}$, lat = ${:.3f}$"
    #           .format(L/1e+3, L/1e+3, lng, lat))
    plt.xlabel('Eastings (m)')
    plt.ylabel('Northings (m)')
    if show:
        plt.show()

def fieldobs_lookup(readings):

    from gascoyne_config import config_layers
    readstr = [ ]
    for i, v in enumerate(readings):
        if v in config_layers.index:
            readstr.append(config_layers.loc[v,'name'])
        else:
            readstr.append('Unknown layer')
    return readstr

def main_contours():
    """
    The main routine to run a suite of diagnostic plots
    """

    # Load everything
    magSensors = np.loadtxt("magSensors.csv", delimiter=',')
    magReadings = np.loadtxt("magReadings.csv", delimiter=',')
    gravSensors = np.loadtxt("gravSensors.csv", delimiter=',')
    gravReadings = np.loadtxt("gravReadings.csv", delimiter=',')
    samples = np.load(pickaxe_npz)
    N = len(samples['magReadings'])

    # Make a few plots of sensors
    plot_sensor(magSensors, magReadings, samples['magReadings'][N/2:],
                units='nT', show=False) #, energy=samples['energy'][N/2:])
    plt.savefig('mag_contours.png')
    plt.show()
    plot_sensor(gravSensors, gravReadings, samples['gravReadings'][N/2:],
                units='mgal', show=False) #, energy=samples['energy'][N/2:])
    plt.savefig('grav_contours.png')
    plt.show()

def main_fieldobs():
    """
    The main routine to show simulated field observations
    """

    # First show the data we expect
    fieldSensors = pd.read_csv('fieldobsSensors.csv', names=['x','y'], comment='#')
    fieldReadings = pd.read_csv('fieldobsReadings.csv', names=['val'], comment='#')
    fieldLabels = fieldSensors.assign(val=fieldobs_lookup(fieldReadings.val))
    fig = plt.figure(figsize=(6,6))
    display_ground_truth(fieldLabels, show=False)
    plt.title('Field Observations')
    plt.savefig('boundary_data.png')

    # Now show samples
    samples = np.load(pickaxe_npz)
    fig = plt.figure(figsize=(6,6))
    i = len(samples['fieldReadings'])
    readings = samples['fieldReadings'][i-1]
    fieldLabels.val = fieldobs_lookup(readings)
    display_ground_truth(fieldLabels, show=False)
    plt.title('Forward-Modeled Field Observations, '
              'Sample {} from MCMC Chain'.format(i))
    plt.savefig('boundary_fwdmodel_endchain.png')
    plt.close()

def main_boundarymovie():
    """
    Makes a movie of how the boundaries change as the chain samples
    """
    # Load everything
    magSensors = np.loadtxt("magSensors.csv", delimiter=',')
    magReadings = np.loadtxt("magReadings.csv", delimiter=',')
    gravSensors = np.loadtxt("gravSensors.csv", delimiter=',')
    gravReadings = np.loadtxt("gravReadings.csv", delimiter=',')
    samples = np.load(pickaxe_npz)

    # Try fitting a few GP layers
    layer_labels = ['layer{}ctrlPoints'.format(i) for i in range(4)]
    for i in np.arange(0, 2500, 25):
        layer_pars = np.array([samples[ll][i] for ll in layer_labels]).reshape(4,5,5)
        fig = plt.figure()
        gp_predict(gravSensors, layer_pars, ((0.0, 2e+4), (0.0, 2e+4)))
        plt.savefig('boundary_movie_frame{:04d}.png'.format(i))
        plt.close()


if __name__ == "__main__":
    main_contours()
    # main_fieldobs()
