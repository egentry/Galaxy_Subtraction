def main():
    """
    Takes in a drizzled fits image (final weight type "EXP") and creates a Poisson photon-noise map. 


    Things to do later:
        Check about taking into account gain

        In units of COUNTS or COUNTS/SEC?
            How does a sci_drz image report data in terms of counts?
                Normalized to nominal exposure time?

        Read in gain, exposure time (don't leave as hardcoded)

   
        


    Version history
        2013.08.12 - Eric Gentry - Began writing script
        2013.08.13 - Eric Gentry - initial draft completed
    """

    from astropy.io import fits
    import numpy as np
    import os

###

    main_dir = '/home/egentry/Data/HST/calib_attempt3/'

    processed_file_dir = 'Processed_Files/'
    

    if not os.path.exists( main_dir + processed_file_dir):
        os.makedirs(main_dir + processed_file_dir)

    sci_filename = 'f160w_drz_sci.fits'
    wht_filename = 'f160w_drz_wht.fits'
    ctx_filename = 'f160w_drz_ctx.fits'


###

    exp_time = 2708.81 # seconds

    sci_hdulist = fits.open(main_dir + sci_filename)
    wht_hdulist = fits.open(main_dir + wht_filename)
    ctx_hdulist = fits.open(main_dir + ctx_filename)

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

    error_data      = np.empty( sci_data.shape )
    bad_pixel_data  = np.zeros( sci_data.shape )

    for x_i in xrange(sci_size_x):
        for y_i in xrange(sci_size_y):
            if ctx_data[y_i][x_i] == 0:
                bad_pixel_data[y_i][x_i] = 1

            if np.isnan(sci_data[y_i][x_i]):
                sci_data[y_i][x_i] = 0
                bad_pixel_data[y_i][x_i] = 1
            
            if np.isinf(sci_data[y_i][x_i]):
                bad_pixel_data[y_i][x_i] = 1
                sci_data[y_i][x_i] = 0

            if wht_data[y_i][x_i] == 0:
                bad_pixel_data[y_i][x_i] = 1 
                wht_data[y_i][x_i] = 1       
    

    error_data = sci_data * wht_data

    error_data = np.absolute(error_data) ** .5
    
    error_data = np.divide(error_data, wht_data)  # elementwise



    bad_pixel_hdulist = fits.PrimaryHDU(bad_pixel_data)
    error_hdulist     = fits.PrimaryHDU(error_data)


    bad_pixel_hdulist.writeto(main_dir + processed_file_dir + 'bad_pixel.fits', clobber=True)
    error_hdulist.writeto(main_dir + processed_file_dir + '/error.fits', clobber=True)


    
main()
