#!/usr/bin/env python

"""
RS 2017/06/22:  Nature Conservancy Urban Green Space + Respiratory Health
    Web Mercator Projection Mapping

This set of functions is supposed to help us query the ArcGIS image servers
for the NSW Aerial Imaging dataset on Australia's NationalMap.  No doubt it
re-invents some wheels, but it's in Python and we finally understand it.
"""

import numpy as np

# ArcGIS coords run to +/- 20,000 km, or about 1/2 Earth's circumference
D = 6378137.0 * np.pi

# Calibration for apparent offset in projection based on Google Earth;
# this is worth querying further, since we don't know where it comes from,
# but it should be good enough for Sydney/NSW in particular.
# The sizes of these coefficients give a rough order of magnitude for
# the overall accuracy of the projection.
x_offset, y_offset = +6.0, -0.5

# Fractional tolerance of projection accuracy in distance
fractol = 1e-10

# exception class
class WebMercatorException(Exception):

    def __init__(self, msg):
        self.msg = msg

    def __str__(self, msg):
        return self.msg

def mercator(lng, lat):
    """
    Calculates (x,y) coordinates in the "Web Mercator" projection currently
    used for the ArcGIS image servers; see
        https://en.wikipedia.org/wiki/Mercator_projection#The_spherical_model
        https://en.wikipedia.org/wiki/Web_Mercator
    Tests using Google Earth lat/lng suggest this is accurate to about 6 m.
        lng, lat:  longitude and latitude of image centre in decimal degrees
    """
    x = D/np.pi * np.radians(lng) + x_offset
    y = D/np.pi * np.log(np.tan(np.pi/4 + np.radians(lat)/2)) + y_offset
    return x, y

def inv_mercator(x, y):
    """
    Calculates (x,y) coordinates in the "Web Mercator" projection currently
    used for the ArcGIS image servers; see
        https://en.wikipedia.org/wiki/Mercator_projection#The_spherical_model
        https://en.wikipedia.org/wiki/Web_Mercator
    Tests using Google Earth lat/lng suggest this is accurate to about 6 m.
        x, y:  Web Mercator projected coordinates in esriMeters
    """
    lng = np.degrees((x - x_offset)*np.pi/D)
    lat = np.degrees(2*np.arctan(np.exp((y - y_offset)*np.pi/D)) - np.pi/2)
    return lng, lat

def bounding_box(L, x=None, y=None, lng=None, lat=None):
    """
    Calculates the corners of a square ArcGIS bounding box.
        L:  side of box, in esriMeters
    """
    if lng and lat:
        x, y = mercator(lng, lat)
    elif x and y:
        pass
    else:
        raise WebMercatorException(
                "Need to specify either lng and lat, or x and y")
    return x - 0.5*L, y - 0.5*L, x + 0.5*L, y + 0.5*L

def test_inv_mercator_1():
    """
    Runs a quick test on the invertibility of our projections.
    Lays down a grid in (lat, lng), projects, deprojects, and compares with
    original grid.  Raises an exception if max deviation exceeds fractol.
    """
    latlist = np.arange(-75, 75.1, 10)
    lnglist = np.arange(0, 360, 10)
    lng, lat = np.meshgrid(lnglist, latlist)
    resids = np.array([lng, lat]) - inv_mercator(*mercator(lng, lat))
    maxfracdev = np.max(np.abs(resids))/180.0
    if maxfracdev > fractol:
        raise WebMercatorException(
                "Max fractional deviation in projection exceeds tolerance")
    return True

def test_inv_mercator_2():
    """
    Runs another quick test on the invertibility of our projections.
    Lays down a grid in (x, y), deprojects, projects, and compares with
    original grid.  Raises an exception if max deviation exceeds fractol.
    """
    xlist = np.linspace(-D, D, 15)
    ylist = np.linspace(-0.9*D, 0.9*D, 15)
    x, y = np.meshgrid(xlist, ylist)
    resids = np.array([x, y]) - mercator(*inv_mercator(x, y))
    maxfracdev = np.max(np.abs(resids))/D
    if maxfracdev > fractol:
        raise WebMercatorException(
                "Max fractional deviation in projection exceeds tolerance")
    return True
