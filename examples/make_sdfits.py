
import numpy as np

from gbt_synth import sample_sky 
from gbt_synth.synth_sdfits import make_synth


if __name__ == '__main__':

    ## Use these cards to create a WCS transformation.
    ## Altough it uses units for the spectral axis, it will not work if you use
    ## something different than velocity in km/s.
    #cards = {'CRVAL1': 80,         'CRVAL2': 0,          'CRVAL3': -100, 
    #         'CTYPE1': 'GLON-SFL', 'CTYPE2': 'GLAT-SFL', 'CTYPE3': 'VELO',
    #         'CDELT1': -5/60.,     'CDELT2': 5/60.,      'CDELT3': 10,
    #         'CRPIX1': 20,         'CRPIX2': 20,         'CRPIX3': 1,
    #         'CUNIT1': 'deg',      'CUNIT2': 'deg',      'CUNIT3': 'km/s',
    #         'OBJECT': 'Cygnus-X'}
    #shape = (20, 100, 100) # z, y, x - to respect how fits cubes are arranged.
    ## Sources is a list of dictionaries.
    ## x0 and y0 are the location of the source in degrees, radi the source size in degrees (it assume a circle), 
    ## type {line, continuum} specifies if to fill all the channels or if it should look like a Gaussian line.
    ## If type is line, then you must also provide z0, central velocity (km/s), and the line fwhm (km/s).
    ## tb gives the brightness of the continuum or line peak in the brightness units you're using.
    #sources = [{'name': 'DR21', 'x0': 81.6807, 'y0': 0.5374, 'radi': 1., 'tb': 0.5, 'z0': 0, 'fwhm':20, 'type': 'line'},
    #           {'name': 'Gamma Cygni SNR', 'x0': 78.2, 'y0': 2.1, 'radi': 1., 'tb': 1, 'type': 'continuum'},
    #          ]
    #noise = 0.1 # Add Gaussian noise with this standard deviation to the data.
    ## nu0 is the rest frequency for the data and must be in Hz.
    #make_synth('cygx', shape, cards, sources, noise, nu0=1420.4e6, save_cube=True) 

    ## Now lets make something closer to the Galactic pole.
    #cards = {'CRVAL1': 80,         'CRVAL2': 88,         'CRVAL3': -100,
    #         'CTYPE1': 'GLON-TAN', 'CTYPE2': 'GLAT-TAN', 'CTYPE3': 'VELO',
    #         'CDELT1': -5/60.,     'CDELT2': 5/60.,      'CDELT3': 10,
    #         'CRPIX1': 25,         'CRPIX2': 25,         'CRPIX3': 1,
    #         'CUNIT1': 'deg',      'CUNIT2': 'deg',      'CUNIT3': 'km/s',
    #         'OBJECT': 'Pole'}
    #shape = (20, 150, 50) 
    #sources = [{'name': '1', 'x0': 81.6807, 'y0': 89.0, 'radi': 1., 'tb': 0.5, 'z0': 0, 'fwhm':20, 'type': 'line'},
    #           {'name': '2', 'x0': 80.0,    'y0': 88.0, 'radi': 1., 'tb': 0.6, 'z0': 0, 'fwhm':30, 'type': 'line'},
    #          ]
    #noise = 0.05 
    #make_synth('gpole', shape, cards, sources, noise, nu0=1420.4e6, save_cube=True)

    #cards = {'CRVAL1': 80,         'CRVAL2': 0,          'CRVAL3': -100, 
    #         'CTYPE1': 'GLON-SFL', 'CTYPE2': 'GLAT-SFL', 'CTYPE3': 'VELO',
    #         'CDELT1': -5/60.,     'CDELT2': 5/60.,      'CDELT3': 10,
    #         'CRPIX1': 20,         'CRPIX2': 20,         'CRPIX3': 1,
    #         'CUNIT1': 'deg',      'CUNIT2': 'deg',      'CUNIT3': 'km/s',
    #         'OBJECT': 'Cygnus-X'}
    #shape = (10, 1000, 10) # z, y, x - to respect how fits cubes are arranged.
    ## Sources is a list of dictionaries.
    ## x0 and y0 are the location of the source in degrees, radi the source size in degrees (it assume a circle), 
    ## type {line, continuum} specifies if to fill all the channels or if it should look like a Gaussian line.
    ## If type is line, then you must also provide z0, central velocity (km/s), and the line fwhm (km/s).
    ## tb gives the brightness of the continuum or line peak in the brightness units you're using.
    #sources = [{'name': 'DR21', 'x0': 81.6807, 'y0': 0.5374, 'radi': 1., 'tb': 0.5, 'z0': 0, 'fwhm':20, 'type': 'line'},
    #           {'name': 'Gamma Cygni SNR', 'x0': 78.2, 'y0': 2.1, 'radi': 1., 'tb': 1, 'type': 'continuum'},
    #          ]
    #noise = 0.1 # Add Gaussian noise with this standard deviation to the data.
    ## nu0 is the rest frequency for the data and must be in Hz.
    #make_synth('tall', shape, cards, sources, noise, nu0=1420.4e6, save_cube=True) 

    #cards = {'CRVAL1': 0,         'CRVAL2': 0,          'CRVAL3': -100, 
    #         'CTYPE1': 'GLON-SFL', 'CTYPE2': 'GLAT-SFL', 'CTYPE3': 'VELO',
    #         'CDELT1': -1/60.,     'CDELT2': 1/60.,      'CDELT3': 10,
    #         'CRPIX1': 50,         'CRPIX2': 50,         'CRPIX3': 1,
    #         'CUNIT1': 'deg',      'CUNIT2': 'deg',      'CUNIT3': 'km/s',
    #         'OBJECT': 'Cygnus-X'}
    #shape = (10, 100, 100) # z, y, x - to respect how fits cubes are arranged.
    #sources = [{'name': 'DR21', 'x0': 0, 'y0': 0, 'radi': 0.0001, 'tb': 10., 'type': 'continuum'},
    #          ]
    #noise = 0.001 # Add Gaussian noise with this standard deviation to the data.
    #make_synth('point', shape, cards, sources, noise, nu0=1420.4e6, save_cube=True) 

    #cards = {'CRVAL1': 0,         'CRVAL2': 0,          'CRVAL3': -100,
    #         'CTYPE1': 'GLON-SFL', 'CTYPE2': 'GLAT-SFL', 'CTYPE3': 'VELO',
    #         'CDELT1': -2/60.,     'CDELT2': 2/60.,      'CDELT3': 10,
    #         'CRPIX1': 50,         'CRPIX2': 50,         'CRPIX3': 1,
    #         'CUNIT1': 'deg',      'CUNIT2': 'deg',      'CUNIT3': 'km/s',
    #         'OBJECT': 'Cygnus-X'}
    #shape = (10, 100, 100) # z, y, x - to respect how fits cubes are arranged.
    #sources = [{'name': 'DR21', 'x0': 0, 'y0': 0, 'radi': 0.0001, 'tb': 10., 'type': 'continuum'},
    #          ]
    #noise = 0.1 # Add Gaussian noise with this standard deviation to the data.
    #make_synth('point_5amin', shape, cards, sources, noise, nu0=1420.4e6, save_cube=True)

    #shape = (10,100,100)
    #cards = {'CRVAL1': 0,          'CRVAL2': 0,          'CRVAL3': -100,
    #         'CTYPE1': 'GLON-TAN', 'CTYPE2': 'GLAT-TAN', 'CTYPE3': 'VELO',
    #         'CDELT1': -2/60.,     'CDELT2': 2/60.,      'CDELT3': 10,
    #         'CRPIX1': shape[2]//2,'CRPIX2': shape[1]//2,'CRPIX3': 1,
    #         'CUNIT1': 'deg',      'CUNIT2': 'deg',      'CUNIT3': 'km/s',
    #         'OBJECT': 'Cygnus-X'}
    ##shape = (10, 100, 100) # z, y, x - to respect how fits cubes are arranged.
    #sources = [{'name': 'W75N', 'x0': 0.5, 'y0': 0.5, 'radi': 0.05, 'tb': 10., 'type': 'continuum'},
    #           {'name': 'DR21', 'x0': 0.5, 'y0': -0.5, 'radi': 0.1, 'tb': 10., 'type': 'continuum'},
    #           {'name': 'line', 'x0': -0.5,'y0': -0.5, 'radi': 0.1, 'tb': 10., 'type': 'line', 'z0': -5, 'fwhm': 5}
    #          ]
    #noise = 0.1 # Add Gaussian noise with this standard deviation to the data.
    #make_synth('tanproj', shape, cards, sources, noise, nu0=1420.4e6, save_cube=True)

#    shape = (2,100,100)
#    restfreq = 7e9
#    hpbw = 0.029446016
#    cards = {'CRVAL1': 10,         'CRVAL2': 10,         'CRVAL3': 0,
#             'CTYPE1': 'GLON-SFL', 'CTYPE2': 'GLAT-SFL', 'CTYPE3': 'VELO',
#             'CDELT1': -0.1/60.,   'CDELT2': 0.1/60.,    'CDELT3': 10,
#             'CRPIX1': shape[2]//2,'CRPIX2': shape[1]//2,'CRPIX3': 1,
#             'CUNIT1': 'deg',      'CUNIT2': 'deg',      'CUNIT3': 'km/s',
#             'NAXIS1': shape[2],   'NAXIS2': shape[1],   'NAXIS3': shape[0],
#             'OBJECT': 'SgrA*',    'RESTFRQ': restfreq,  
#             'BMAJ'  : hpbw,       'BMIN'   : hpbw,      'BPA': 0}
#
#    sources = [{'name': 'SgrA*', 'x0': 10, 'y0': 10, 'radi': 0.0005, 'tb': 10., 'type': 'continuum'}]
#    noise = 0.001
#    make_synth('sgra', shape, cards, sources, noise, nu0=restfreq, save_cube=True, convolve_cube=True, hpbw=hpbw)
#

    # Galaxy with a HI disk.
    shape = (200,250,250)
    restfreq = 1.42e9
    lmbd = 3e8/restfreq
    hpbw = np.rad2deg(1.2*lmbd/100.)
    loc = [140.51115792, 50.97659364]
    cards = {'CRVAL1': loc[0],     'CRVAL2': loc[1],     'CRVAL3': 0,
             'CTYPE1': 'RA---SFL', 'CTYPE2': 'DEC--SFL', 'CTYPE3': 'VELO',
             'CDELT1': -0.2/60.,   'CDELT2': 0.2/60.,    'CDELT3': 10,
             'CRPIX1': shape[2]//2,'CRPIX2': shape[1]//2,'CRPIX3': 100,
             'CUNIT1': 'deg',      'CUNIT2': 'deg',      'CUNIT3': 'km/s',
             'NAXIS1': shape[2],   'NAXIS2': shape[1],   'NAXIS3': shape[0],
             'OBJECT': 'NGC2841',  'RESTFRQ': restfreq,  
             'BMAJ'  : hpbw,       'BMIN'   : hpbw,      'BPA': 0}
    sources = [#{'name': 'NGC2841', 'x0': loc[0], 'y0': loc[1], 'radi': 3/60., 'tb': 10., 'type': 'continuum'},
               {'name': 'HI-disk', 'x0': loc[0], 'y0': loc[1], 'radi': 3/60., 'tb': 0.5, 'type': 'line', 'z0': 0, 'fwhm': 220},
               {'name': 'HI-edisk', 'x0': loc[0], 'y0': loc[1], 'radi': 20/60., 'tb': 0.2, 'type': 'line', 'z0': 100, 'fwhm': 220},
               {'name': 'HI-edisk', 'x0': loc[0], 'y0': loc[1], 'radi': 20/60., 'tb': 0.1, 'type': 'line', 'z0': -90, 'fwhm': 220},
              ]
    noise = 0.01
    make_synth('NGC2841', shape, cards, sources, noise, nu0=restfreq, save_cube=True, convolve_cube=True, hpbw=hpbw)
    sample_sky.main("NGC2841_tracks_sdfits.fits", "NGC2841_cube.fits", 
                    "/home/scratch/psalas/user_support/AGBT22A_287/cats/ngc2841.csv",
                    cards, nu0=restfreq)
