
def main():
    """
    Takes in an existing bad pixel mask, and additionally masks out inner radii of galaxy

   
        


    Version history
        2013.10.17 - Eric Gentry - initial draft completed
    """

    from astropy.io import fits
    import numpy as np

###

    main_dir = '/home/egentry/Data/HST/PictorA/bad_pixel_masks/'

    mask_filename_in  = 'f160w_drz_mask.fits'
    mask_filename_out = 'f160w_drz_mask_OuterRadii.fits'


###

    mask_in_hdulist = fits.open(main_dir + mask_filename_in)

    mask_in_data = mask_in_hdulist[0].data


    mask_size_x = mask_in_data.shape[1]
    mask_size_y = mask_in_data.shape[0]

    #Center region to crop out:
    center_radius = 75  # pixels
    center_x      = 758 #pixel, zero-indexed
    center_y      = 265 #pixel, zero-indexed

    mask_out_data       = np.zeros( mask_in_data.shape )

    for x_i in xrange(mask_size_x):
        for y_i in xrange(mask_size_y):
            if (center_x - x_i)**2 + (center_y - y_i)**2 < center_radius**2:
                mask_out_data[y_i][x_i] = 1
            else:
                mask_out_data[y_i][x_i] = mask_in_data[y_i][x_i]
                 



    mask_out_hdulist = fits.PrimaryHDU(mask_out_data)

    mask_out_hdulist.writeto(main_dir + mask_filename_out, clobber=True)


    

