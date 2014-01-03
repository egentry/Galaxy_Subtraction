def create_bad_pixel_mask(input_dir, output_dir, filter_name):
    """
    TO DO: incorporate read noise error estimate.
            Currently GALFIT does a better job at estimating error -- this is only useful for making bad pixel map   
        


    Version history
        2013.08.12 - Eric Gentry - Began writing script
        2013.08.13 - Eric Gentry - initial draft completed
        2014.01.02 - Eric Gentry - Removed 'Error Image' feature - GALFIT -outsig was better
    """
    from astropy.io import fits
    import numpy as np
    import os

###
    

    if not os.path.exists( output_dir):
        os.makedirs(output_dir)

    sci_filename = filter_name + '_drz_sci.fits'
    wht_filename = filter_name + '_drz_wht.fits'
    ctx_filename = filter_name + '_drz_ctx.fits'


###


    sci_hdulist = fits.open(input_dir + sci_filename)
    wht_hdulist = fits.open(input_dir + wht_filename)
    ctx_hdulist = fits.open(input_dir + ctx_filename)

    sci_data = sci_hdulist[0].data
    wht_data = wht_hdulist[0].data
    ctx_data = ctx_hdulist[0].data

    sci_size_x = sci_data.shape[1]
    sci_size_y = sci_data.shape[0]

    wht_size_x = wht_data.shape[1]
    wht_size_y = wht_data.shape[0]

    ctx_size_x = ctx_data.shape[1]
    ctx_size_y = ctx_data.shape[0]

    if ( sci_size_x != wht_size_x ) or ( sci_size_y != wht_size_y):
        print 'SCI and WHT dimensions do not match'
        print 'SCI Dimensions: ', sci_size_x, ' by ' , sci_size_y
        print 'SCI Dimensions: ', wht_size_x, ' by ' , wht_size_y
        return

    if ( sci_size_x != ctx_size_x ) or ( sci_size_y != ctx_size_y):
        print 'SCI and CTX dimensions do not match'
        print 'SCI Dimensions: ', sci_size_x, ' by ' , sci_size_y
        print 'CTX Dimensions: ', ctx_size_x, ' by ' , ctx_size_y
        return

    # bad_pixel_data  = np.zeros( sci_data.shape )

    bad_pixel_data = np.logical_not( np.isfinite(sci_data) * ctx_data * wht_data) * 1

    # # set infinite pixels to zero, and set corresponding wht, ctx values to 0
    # ctx_data = ctx_data * np.isfinite(sci_data)
    # wht_data = wht_data * np.isfinite(sci_data)
    # sci_data = sci_data * np.isfinite(sci_data)




    bad_pixel_hdulist = fits.PrimaryHDU(bad_pixel_data)
    bad_pixel_hdulist.writeto(output_dir + filter_name + '_bad_pixel.fits', clobber=True)

    # sci_hdulist.writeto(input_dir + sci_filename, clobber=True)
    # wht_hdulist.writeto(input_dir + wht_filename, clobber=True)
    # ctx_hdulist.writeto(input_dir + ctx_filename, clobber=True)

    print filter_name, ' pixel mask done \n'

    
def main_f160w():

    input_dir   = '/home/egentry/Data/HST/PictorA/IBJX01010/drizzled/'
    output_dir  = '/home/egentry/Data/HST/PictorA/bad_pixel_masks/'
    filter_name = 'f160w'

    create_bad_pixel_mask(input_dir, output_dir, filter_name)

def main_f475w():
    input_dir   = '/home/egentry/Data/HST/PictorA/IBJX01020/drizzled/'
    output_dir  = '/home/egentry/Data/HST/PictorA/bad_pixel_masks/'
    filter_name = 'f475w'

    create_bad_pixel_mask(input_dir, output_dir, filter_name)

def main_f814w():
    input_dir   = '/home/egentry/Data/HST/PictorA/IBJX01030/drizzled/'
    output_dir  = '/home/egentry/Data/HST/PictorA/bad_pixel_masks/'
    filter_name = 'f814w'

    create_bad_pixel_mask(input_dir, output_dir, filter_name)

main_f160w()
main_f475w()
main_f814w()


    






