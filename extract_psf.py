def extract(drz_dir, filter_name, psf_center, psf_width, psf_height):
    """ Extracts a PSF from a particular star in a particular dataset

        Expects a DRZ file from the HST Proposal 12261, sampled at 2x the resolution of the FLT 
        
        Pixels here follow a zero indexing convention. Be careful -- DS9 uses 1 indexing!

        Currently this is only meant for a very specific purpose. Many values are hardcoded. This may be fixed later.


        Version History
            2013.08.30 -- Eric Gentry -- Created limited purpose draft

    """

    from astropy.io import fits
    import numpy as np
    import os

    #0-indexed:


    min_psf_x = psf_center[0] - (psf_width / 2)
    max_psf_x = psf_center[0] + (psf_width / 2)
    min_psf_y = psf_center[1] - (psf_height / 2)
    max_psf_y = psf_center[1] + (psf_height / 2)

    # print min_psf_x
    # print min_psf_y
    # print max_psf_x
    # print max_psf_y

    
    drz_filename = drz_dir + filter_name + '_drz_sci.fits'
    psf_dir      = '/home/egentry/Data/HST/PictorA/PSFs/' + filter_name + '/extracted/' # to be created
    # psf_medsub_filename = 'psf_medsub_extracted_test.fits' # to be created

    if not os.path.exists(psf_dir):
        os.makedirs(psf_dir)

    
    drz_hdulist = fits.open(drz_filename)

    drz_data = drz_hdulist[0].data

    psf_data = drz_data[min_psf_y : max_psf_y + 1, min_psf_x : max_psf_x + 1]

    psf_median = np.median(psf_data)

    psf_data = psf_data - psf_median
    psf_data = psf_data.clip(min=0)


    psf_data_hdu = fits.PrimaryHDU(psf_data)
    psf_data_hdu.writeto(psf_dir + filter_name + '_psf_drz_sci_extracted.fits', clobber=True)

    
    drz_hdulist.close()



def main_f160w():

    drz_dir = '/home/egentry/Data/HST/PictorA/IBJX01010/drizzled/'
    filter_name = 'f160w'

    psf_center = (277, 483)

    psf_width  = 288
    psf_height = 244


    extract(drz_dir, filter_name, psf_center, psf_width, psf_height)



def main():
    main_f160w()
    # main_f475w()
    # main_f814w()  

main()
   