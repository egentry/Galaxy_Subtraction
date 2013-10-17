main_dir = '/home/egentry/Data/HST/PictorA/PSFs/extract/'


from astropy.io import fits
import numpy as np
import os


def e():
    """ Extracts a PSF from a particular star in a particular dataset

        Expects a DRZ file from the HST Proposal 12261, sampled at 2x the resolution of the FLT 
        
        Pixels here follow a zero indexing convention. Be careful -- DS9 uses 1 indexing!

        Currently this is only meant for a very specific purpose. Many values are hardcoded. This may be fixed later.


        Version History
            2013.08.30 -- Eric Gentry -- Created limited purpose draft

    """



    #0-indexed:
    psf_center = (277, 483)

    psf_width  = 288
    psf_height = 244

    min_psf_x = psf_center[0] - (psf_width / 2)
    max_psf_x = psf_center[0] + (psf_width / 2)
    min_psf_y = psf_center[1] - (psf_height / 2)
    max_psf_y = psf_center[1] + (psf_height / 2)

    print min_psf_x
    print min_psf_y
    print max_psf_x
    print max_psf_y
    

    
    drz_filename = 'f160w_drz_sci.fits'
    psf_filename = 'psf_extracted.fits' # to be created
    
    drz_hdulist = fits.open(main_dir + drz_filename)

    drz_data = drz_hdulist[0].data

    psf_data = drz_data[min_psf_y : max_psf_y + 1, min_psf_x : max_psf_x + 1]

    print drz_data
    print psf_data
    print drz_data[psf_center[1], psf_center[0]]

    psf_data_hdu = fits.PrimaryHDU(psf_data)
    psf_data_hdu.writeto(main_dir + psf_filename, clobber=True)
    
    drz_hdulist.close()



    
