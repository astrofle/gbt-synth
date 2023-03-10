"""
Synthetic SDFITS writer.
For now, it only creates circular sources.
"""

import numpy as np

import fitsio
import astropy.constants as ac

from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import convolve_fft, Gaussian2DKernel


from groundhog import sd_fits_keys


def angular_distance(l1, b1, l2, b2):
    """
    """

    return 2. * np.arcsin(np.sqrt(
        np.sin((b1 - b2) / 2.) ** 2 +
        np.cos(b1) * np.cos(b2) * np.sin((l1 - l2) / 2.) ** 2
        ))


def get_axes_wcs(wcs, nx, ny, nz):
    """
    Generate spatial and spectral axes from a WCS object.
    """

    z = wcs.spectral.all_pix2world(np.arange(nz), 0)[0]*1e-3 # For some reason the WCS object changes to m/s.
    yy,xx = np.mgrid[0:ny:1,0:nx:1]
    pixls = np.vstack((xx.flatten(), yy.flatten())).T
    x, y = wcs.celestial.all_pix2world(xx.flatten(), yy.flatten(), 0)

    return x, y, z


def make_cube(x, y, z, nx, ny, nz, sources, noise):
    """
    """

    cube = np.zeros((nz,ny,nx), dtype=np.float32)

    # Add some noise, with 0.1 K rms
    cube += np.random.normal(loc=0.0, scale=noise, size=cube.shape)

    for source in sources:
        mask, tant = make_source(source, x, y, z)
        cube[:,mask.reshape(ny,nx)] += tant[:,np.newaxis]

    return cube


def make_wcs_header(cards):
    """
    cards : dict
        Dictionary with WCS header keywords.
        E.g., CDELT1, CRVAL1, etc...
    """

    # Create a header and use it to generate our coordinates.
    hdu = fits.PrimaryHDU()
    head = hdu.header
    for k in cards:
        head.set(k, cards[k])

    # Create a WCS object to handle coordinates.
    wcs = WCS(head)

    return wcs


def make_sdfits(output, data, x, y, z, cards, nu0=1420.4e6):
    """
    Write and SDFITS file to `output` and fill it with `data`
    giving it spatial coordinates `x` and `y` and spectral coordinates `z`.
    """

    nrows = data.shape[0]
    nchan = data.shape[1]

    # Define some SDFITS keywords.
    sdfits_keys = {'TELESCOP': 'NRAO_GBT', # The one and only.
                   'INSTRUME': 'VEGAS', # Everyone's favorite backend.
                   'DATE': '2022-03-25',
                   'PROJID': 'TEST_MATRIX',
                   'SITELONG': -79.83983,
                   'SITELAT': 38.43312,
                   'SITEELEV': 824.595,
                   'TFORM7': f'{nchan}E'
                  }

    # Create the primary HDU for the SDFITS and fill it.
    hlist = []
    for k in sd_fits_keys.phdu:
        # If the keyword has a variable value, set it to what we defined above.
        if k in sd_fits_keys.phdu_fill:
            value = sdfits_keys[k]
        else:
            value = sd_fits_keys.phdu[k][0]
        hlist.append({'name':k, 'value':value, 'comment':sd_fits_keys.phdu[k][1]})

    outfits = fitsio.FITS(output, 'rw', clobber=True)
    outfits.write(None, header=hlist)
        
    # Now lets create a table with the data.
    dtype_list = []
    for i,k in enumerate(sd_fits_keys.table_dtypes):
        if k[0] != 'DATA':
            dtype_list.append(k)
        else:
            dtype_list.append((k[0], k[1], nchan,))
            
    dtype = np.dtype(dtype_list)
    new_table = np.zeros(nrows, dtype=dtype)
    # Fill the new table.
    new_table['DATA'] = data
    new_table['OBJECT'] = cards['OBJECT']
    # The gridder only knows about frequencies.
    new_table['CTYPE1'] = 'FREQ-OBS'
    new_table['CRVAL1'] = z[0] # HI
    new_table['CRPIX1'] = 1 # Python starts at 0, FITS at 1
    new_table['CDELT1'] = np.mean(np.diff(z)) # channel width in Hz
    # Position.
    new_table['CTYPE2'] = cards['CTYPE1'].split('-')[0]
    new_table['CTYPE3'] = cards['CTYPE2'].split('-')[0]
    new_table['CRVAL2'] = x
    new_table['CRVAL3'] = y
    # Avoid divisions by zero.
    new_table['TSYS'] = 20.
    new_table['EXPOSURE'] = 1.
    # Other stuff to avoid errors with the gridder.
    new_table['RESTFREQ'] = nu0
    new_table['DATE-OBS'] = '2022-03-25T18:41:09.50'
    new_table['FRONTEND'] = 'Rcvr1_2'
    new_table['VELDEF'] = 'RADI-LSR'

    # Add a table extension to the sdfits.
    outfits.create_table_hdu(dtype=new_table.dtype, extname='SINGLE DISH')

    # Make the table header.
    hlist = []
    for row in sd_fits_keys.hdu:
        # Give the correct number of channels.
        if row[0] == 'TFORM7':
            value = f'{nchan}E'
        else:
            value = row[1]
        hlist.append({'name':row[0], 'value':value, 'comment':row[2]})
    outfits[-1].write_keys(hlist)
    outfits[-1]._update_info()

    # Add the table.
    outfits[-1].append(new_table)

    # Close.
    outfits.close()


def make_source(source, x, y, z):
    """
    """
    
    dist = angular_distance(np.deg2rad(source['x0']), 
                            np.deg2rad(source['y0']), 
                            np.deg2rad(x), 
                            np.deg2rad(y))
    mask = abs(dist) <= np.deg2rad(source['radi'])

    # Add a Gaussian line.
    if source['type'] == 'line':
        sigma = source['fwhm']/(2.*np.sqrt(2.*np.log(2.)))
        tant = source['tb']*np.exp(-(z - source['z0'])**2/(2*sigma**2.))
    elif source['type'] == 'continuum':
        tant = np.ones(len(z))*source['tb']

    return mask, tant


def smooth_cube(cube, hpbw, wcs):
    """
    """

    cube_smo = np.empty_like(cube)

    sigma_y = hpbw/abs(wcs.wcs.cdelt[1])/(np.sqrt(8*np.log(2)))
    sigma_x = hpbw/abs(wcs.wcs.cdelt[0])/(np.sqrt(8*np.log(2)))
    kernel = Gaussian2DKernel(x_stddev=sigma_x, y_stddev=sigma_y)

    for i in range(wcs.array_shape[0]):
        cube_smo[i] = convolve_fft(cube[i], kernel, preserve_nan=True)

    return cube_smo


def make_synth(basename, shape, wcs_dict, sources, noise, nu0=1420.4e6, 
               convolve_cube=False, hpbw=0, save_cube=True):
    """
    """

    # Create a WCS object given the user provided keywords.
    wcs = make_wcs_header(wcs_dict)
    # Get the axes values in WCS for the given shape.
    x, y, z = get_axes_wcs(wcs, shape[2], shape[1], shape[0])
    # Fill the cube with sources and noise.
    cube = make_cube(x, y, z, shape[2], shape[1], shape[0], sources, noise)
    if convolve_cube:
        # Convolve the cube to degrade its spatial resolution.
        cube = smooth_cube(cube, hpbw, wcs)
    if save_cube:
        fits.writeto(f'{basename}_cube.fits', cube, header=wcs.to_header())
    # Now, lets make an sdfits file.
    data = cube.reshape((shape[0],shape[1]*shape[2])).T
    # Convert the spectral axis to frequency.
    z_nu = nu0*(1. - z/ac.c.to(wcs_dict['CUNIT3']).value) # Hz
    output = f'{basename}_sdfits.fits'
    make_sdfits(output, data, x, y, z_nu, wcs_dict, nu0=nu0)

# Start of main program.
#if __name__ == "__main__":
