main_dir = '/home/egentry/Data/HST/PictorA/PSFs/drizzler/'

artificial_star_dir = main_dir + 'artificial_stars/'
drizzled_dir        = main_dir + 'drizzled/'
cropped_dir         = main_dir + 'cropped/'
psf_rootname        = 'psf'


from astropy.io import fits
import os


def prepare():
    """ Prepare PSFs to be drizzled 

        Expects FLT + PSF fits files
            Number of FLT files should match number of PSF files
            FLT files should be multiextension .fits files, as given by NASA ('SCI' = index 1)

        Pixels here follow a zero indexing convention. Be careful -- DS9 uses 1 indexing!

        Currently this is only meant for a very specific purpose. Many values are hardcoded. This may be fixed later.


        Version History
            2013.08.09 -- Eric Gentry -- Created limited purpose draft

    """

    if not os.path.exists(artificial_star_dir):
        os.makedirs(artificial_star_dir)



    #0-indexed:
    agn_centers = [ (363, 143), (366, 146) , (369, 150)]

    psf_size = 276
    psf_center = 138


    for Image_index in xrange(3):

        #0 indexed:
        min_nonzero_x = agn_centers[Image_index][0] - psf_center
        max_nonzero_x = agn_centers[Image_index][0] + psf_center -1
        min_nonzero_y = agn_centers[Image_index][1] - psf_center
        max_nonzero_y = agn_centers[Image_index][1] + psf_center -1

    
        psf_filename =         'psf_{}.fits'.format(Image_index)
        flt_filename = 'ibjx01y{}q_flt.fits'.format(Image_index)  # To be done later: accept more general filenames
        
        psf_hdulist = fits.open(main_dir + psf_filename)
        
        flt_hdulist = fits.open(main_dir + flt_filename)

        
        flt_data = flt_hdulist[1].data
        psf_data = psf_hdulist[0].data

        flt_size_x = flt_data.shape[1]
        flt_size_y = flt_data.shape[0]

        psf_size_x = psf_data.shape[1]
        psf_size_y = psf_data.shape[0]

        
        flt_data.fill(0)

#        for x_i in xrange(flt_size_x):
#            for y_i in xrange(flt_size_y):
#                flt_data[y_i][x_i] =  0

        for x_i in xrange(psf_size_x):
            for y_i in xrange(psf_size_y):
                flt_data[y_i + min_nonzero_y][x_i + min_nonzero_x] = psf_data[y_i][x_i] * 1e6

        flt_hdulist.writeto(artificial_star_dir + flt_filename, clobber=True)
        flt_hdulist.close()
        psf_hdulist.close()

def drizzle(clean_files=True):
    """ 

    """
    from drizzlepac import astrodrizzle
    import glob

    if not os.path.exists(drizzled_dir):
        os.makedirs(drizzled_dir)

    

    astrodrizzle.AstroDrizzle(artificial_star_dir + '*flt.fits', output=drizzled_dir + 'psf', static=False, skysub=False,
        driz_separate=False, median=False, blot=False, driz_cr=False, 
        driz_combine=True, final_wht_type='EXP', final_pixfrac=0.7, final_bits=4096,
        final_wcs=True, final_scale=0.06666)

    if clean_files == True:
        for filename in glob.glob(artificial_star_dir +'*_final_mask.fits'):
            os.remove(filename)
                        


def crop():
    """
        Takes a drizzled TinyTim PSF, filters out NaN values, crops to 30arcsec box suitable for GALFIT
    """
    
    import numpy as np

    if not os.path.exists(cropped_dir):
        os.makedirs(cropped_dir)



    psf_filename = psf_rootname + '_drz_sci'

    psf_hdulist = fits.open(drizzled_dir + psf_filename + '.fits')

    psf_data = psf_hdulist[0].data

    psf_size_x = psf_data.shape[1]
    psf_size_y = psf_data.shape[0]

    
    # Filter out bad data -- later check why this happens (maybe final_fillval?)
#    for x_i in xrange(psf_size_x):
#        for y_i in xrange(psf_size_y):
#            if np.isnan(psf_data[y_i][x_i]):
#                print 'found a nan at:', x_i, y_i

    # Filter out bad data -- later check why this happens (maybe final_fillval?)
    #Set all NaN's to 0 (Also changes +/- Inf, but that shouldn't affect this data)
    psf_data = np.nan_to_num(psf_data)

    valid_coords =  np.nonzero(psf_data)


    min_nonzero_x = np.amin( valid_coords[1] )
    max_nonzero_x = np.amax( valid_coords[1] )
    min_nonzero_y = np.amin( valid_coords[0] )
    max_nonzero_y = np.amax( valid_coords[0] )

    psf_data_cropped = psf_data[min_nonzero_y : max_nonzero_y + 1, min_nonzero_x : max_nonzero_x + 1] * 1e-6

    psf_cropped_hdulist = fits.PrimaryHDU(psf_data_cropped)
    psf_cropped_hdulist.writeto(cropped_dir + psf_filename + '_cropped.fits', clobber=True)
    
    psf_hdulist.close()


def main():

    prepare()
    drizzle()
    crop()
   

    
