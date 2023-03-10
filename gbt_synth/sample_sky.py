"""
Tools for sampling a sky model at given location.
"""

import numpy as np
import pandas as pd

from spectral_cube import SpectralCube
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy import constants as ac
from astropy.coordinates import SkyCoord

from .synth_sdfits import make_sdfits


def coo2ds9(coo, radius):
    """
    """

    coostr = ",".join(coo.to_string().split(" "))
    frame = coo.frame.name
    ds9str = f"""{frame};circle({coostr},{radius}")"""

    return ds9str


def locs2coords(locs):
    """
    """

    coords = SkyCoord(locs["RA"], locs["DEC"], 
                      unit=(u.hourangle, u.deg), frame="fk5")

    return coords


def read_sky(filename):
    """
    """

    hdu = fits.open(filename)
    data = hdu[0].data
    head = hdu[0].header
    wcs = WCS(head)

    sky = {"data": data,
           "head": head,
           "wcs" : wcs,
          }

    return sky


def read_locs(filename):
    """
    Example file with locations:
    NAME,RA,DEC,RADIUS
    NGC2841,9:22:02.6779,+50:58:35.7371,300
    """

    df = pd.read_csv(filename, delimiter=",", skiprows=0)

    return df


def sample_sky(skymod, locs):
    """
    """

    cube = SpectralCube(skymod["data"], skymod["wcs"])
    
    sddata = np.empty((len(locs), skymod["head"]["NAXIS3"]), dtype=float)
    sdcoos = locs2coords(locs)

    for i,coo in enumerate(sdcoos):

        ds9str = coo2ds9(coo, locs["RADIUS"][i])
        sddata[i] = cube.subcube_from_ds9region(ds9str).mean(axis=(1,2))

    specax = cube.spectral_axis

    return sdcoos.ra.value, sdcoos.dec.value, specax, sddata


def main(output, sky_filename, locs_filename, cards, nu0=1420.4e6):
    """
    """

    skymod = read_sky(sky_filename)
    locs = read_locs(locs_filename)
    x, y, z, data = sample_sky(skymod, locs)
    z = z.to(cards["CUNIT3"]).value
    z_nu = nu0*(1. - z/ac.c.to(cards["CUNIT3"]).value) # Hz

    make_sdfits(output, data, x, y, z_nu, cards, nu0=nu0)
