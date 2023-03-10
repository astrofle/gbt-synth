
import numpy as np

from gbt_synth import sample_sky 
from gbt_synth.synth_sdfits import make_synth


if __name__ == '__main__':

    # Galaxy with a HI disk.
    shape = (200,250,250)              # z, y, x - to respect how fits cubes are arranged.
    restfreq = 1.42e9                  # Rest frequency of the line of interest in Hz.
    lmbd = 3e8/restfreq                # Wavelength in m.
    hpbw = np.rad2deg(1.2*lmbd/100.)   # Beam size in degrees.
    loc = [140.51115792, 50.97659364]  # Location of the simulated galaxy.
    # Use these cards to create a WCS transformation. This will be used to generate a sky model.
    # Altough it uses units for the spectral axis, it will not work if you use
    # something different than velocity in km/s.
    cards = {'CRVAL1': loc[0],     'CRVAL2' : loc[1],     'CRVAL3': 0,
             'CTYPE1': 'RA---SFL', 'CTYPE2' : 'DEC--SFL', 'CTYPE3': 'VELO',
             'CDELT1': -0.2/60.,   'CDELT2' : 0.2/60.,    'CDELT3': 10,
             'CRPIX1': shape[2]//2,'CRPIX2' : shape[1]//2,'CRPIX3': 100,
             'CUNIT1': 'deg',      'CUNIT2' : 'deg',      'CUNIT3': 'km/s',
             'NAXIS1': shape[2],   'NAXIS2' : shape[1],   'NAXIS3': shape[0],
             'OBJECT': 'NGC2841',  'RESTFRQ': restfreq,  
             'BMAJ'  : hpbw,       'BMIN'   : hpbw,       'BPA': 0}
    # Sources is a list of dictionaries.
    # x0 and y0 are the location of the source in degrees, radi the source size in degrees (it assume a circle), 
    # type {line, continuum} specifies if to fill all the channels or if it should look like a Gaussian line.
    # If type is line, then you must also provide z0, central velocity (km/s), and the line fwhm (km/s).
    # tb gives the brightness of the continuum or line peak in the brightness units you're using.
    sources = [#{'name': 'NGC2841', 'x0': loc[0], 'y0': loc[1], 'radi': 3/60., 'tb': 10., 'type': 'continuum'},
               {'name': 'HI-disk', 'x0': loc[0], 'y0': loc[1], 'radi': 3/60., 'tb': 0.5, 'type': 'line', 'z0': 0, 'fwhm': 220},
               {'name': 'HI-edisk', 'x0': loc[0], 'y0': loc[1], 'radi': 20/60., 'tb': 0.2, 'type': 'line', 'z0': 100, 'fwhm': 220},
               {'name': 'HI-edisk', 'x0': loc[0], 'y0': loc[1], 'radi': 20/60., 'tb': 0.1, 'type': 'line', 'z0': -90, 'fwhm': 220},
              ]
    noise = 0.01 # Add Gaussian noise in K.

    # Create the sky model. It will also create an SDFITS assuming a fully sampled observation.
    make_synth('NGC2841', shape, cards, sources, noise, nu0=restfreq, save_cube=True, convolve_cube=True, hpbw=hpbw)
    # Use the sky model to simulate an observation of the regions defined in ngc2841.csv.
    sample_sky.main("NGC2841_tracks_sdfits.fits", "NGC2841_cube.fits", "ngc2841.csv",
                    cards, nu0=restfreq)
