""" Things to do:
        Streamline prepare(), drizzle(), crop() into one function, which accepts arguments for each filter
"""


from astropy.io import fits
import os
import subprocess
from itertools import izip




def prepare_f160w():
    """ Prepare PSFs to be drizzled 

        Expects FLT + PSF fits files
            Number of FLT files should match number of PSF files
            FLT files should be multiextension .fits files, as given by NASA ('SCI' = index 1)

        Pixels here follow a zero indexing convention. Be careful -- DS9 uses 1 indexing!

        Currently this is only meant for a very specific purpose. Many values are hardcoded. This may be fixed later.


        Version History
            2013.08.09 -- Eric Gentry -- Created limited purpose draft

    """

    main_dir = '/home/egentry/Data/HST/PictorA/PSFs/drizzler/'

    artificial_star_dir = main_dir + 'artificial_stars/'
    drizzled_dir        = main_dir + 'drizzled/'
    cropped_dir         = main_dir + 'cropped/'
    psf_rootname        = 'psf'

    if not os.path.exists(artificial_star_dir):
        os.makedirs(artificial_star_dir)



    #0-indexed: [ (x1,y1), (x2,y2), (x3, y3)]
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

def prepare_f475w():
    """ 
        Adaptation of prepare_f160w
    """

    main_dir = '/home/egentry/Data/HST/PictorA/PSFs/f475w/'  

    flt_dir  = '/home/egentry/Data/HST/PictorA/IBJX01020/'

    artificial_star_dir = main_dir + 'artificial_stars/'
    drizzled_dir        = main_dir + 'drizzled/'
    cropped_dir         = main_dir + 'cropped/'
    psf_rootname        = 'psf'

    if not os.path.exists(artificial_star_dir):
        os.makedirs(artificial_star_dir)



    #0-indexed: [ (x1,y1), (x2,y2), (x3, y3)]
    agn_centers = [ (1698, 837), (1700, 840) , (1703, 841)]
    agn_x = [1699, 1701, 1704]
    agn_y = [838, 841, 842]
    
    flt_indices = [4, 6, 8]

    psf_size = 299
    psf_center = 150


    for (x,y, i) in izip(agn_x, agn_y, flt_indices):

        #0 indexed:
        min_nonzero_x = x-1 - psf_center
        max_nonzero_x = x-1 + psf_center -1
        min_nonzero_y = y-1 - psf_center
        max_nonzero_y = y-1 + psf_center -1

    
        psf_filename =         '{0}_{1}.fits'.format(x,y)
        flt_filename = 'ibjx01y{}q_flt.fits'.format(i)  # To be done later: accept more general filenames
        
        psf_hdulist = fits.open(main_dir + psf_filename)
        
        flt_hdulist = fits.open(flt_dir + flt_filename)

        
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

def prepare_f814w():
    """ 
        Adaptation of prepare_f160w
    """

    main_dir = '/home/egentry/Data/HST/PictorA/PSFs/f814w/'  

    flt_dir  = '/home/egentry/Data/HST/PictorA/IBJX01030/'

    artificial_star_dir = main_dir + 'artificial_stars/'
    drizzled_dir        = main_dir + 'drizzled/'
    cropped_dir         = main_dir + 'cropped/'
    psf_rootname        = 'psf'

    if not os.path.exists(artificial_star_dir):
        os.makedirs(artificial_star_dir)



    agn_x = [1699, 1701, 1704]
    agn_y = [838, 841, 843]
    
    flt_indices = ['a', 'c', 'e']

    psf_size = 299
    psf_center = 150


    for (x,y, i) in izip(agn_x, agn_y, flt_indices):

        #0 indexed:
        min_nonzero_x = x-1 - psf_center
        max_nonzero_x = x-1 + psf_center -1
        min_nonzero_y = y-1 - psf_center
        max_nonzero_y = y-1 + psf_center -1

    
        psf_filename =         '{0}_{1}.fits'.format(x,y)
        flt_filename = 'ibjx01y{}q_flt.fits'.format(i)  # To be done later: accept more general filenames
        
        psf_hdulist = fits.open(main_dir + psf_filename)
        
        flt_hdulist = fits.open(flt_dir + flt_filename)

        
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

def drizzle_f160w(clean_files=True):
    """ 

    """

    main_dir = '/home/egentry/Data/HST/PictorA/PSFs/drizzler/'

    artificial_star_dir = main_dir + 'artificial_stars/'
    drizzled_dir        = main_dir + 'drizzled/'
    cropped_dir         = main_dir + 'cropped/'
    psf_rootname        = 'psf'


    from drizzlepac import astrodrizzle
    import glob

    if not os.path.exists(drizzled_dir):
        os.makedirs(drizzled_dir)

    

    astrodrizzle.AstroDrizzle(artificial_star_dir + '*flt.fits', output=drizzled_dir + 'psf', runfile=drizzled_dir + 'astrodrizzle.log',
        preserve=False, static=False, skysub=True, driz_separate=False, median=False, blot=False, driz_cr=False, 
        driz_combine=True, final_wht_type='EXP', final_pixfrac=0.7, final_bits=4096,
        final_wcs=True, final_scale=0.06666)

    if clean_files == True:
        for filename in glob.glob(artificial_star_dir +'*_final_mask.fits'):
            os.remove(filename)


def drizzle_f475w(clean_files=True):
    """ 

    """
    from drizzlepac import astrodrizzle
    import glob



    main_dir = '/home/egentry/Data/HST/PictorA/PSFs/f475w/'  

    flt_dir  = '/home/egentry/Data/HST/PictorA/IBJX01020/'

    artificial_star_dir = main_dir + 'artificial_stars/'
    drizzled_dir        = main_dir + 'drizzled/'
    psf_rootname        = 'psf'

    if not os.path.exists(drizzled_dir):
        os.makedirs(drizzled_dir)

    

    astrodrizzle.AstroDrizzle(artificial_star_dir + '*flt.fits', output=drizzled_dir + 'psf', runfile=drizzled_dir + 'astrodrizzle.log',
        preserve=False, static=False, skysub=True, driz_separate=False, median=False, blot=False, driz_cr=False, 
        driz_combine=True, final_wht_type='EXP', final_pixfrac=0.7, final_bits=4096,
        final_wcs=True, final_scale=0.01981)

    if clean_files == True:
        for filename in glob.glob(artificial_star_dir +'*_final_mask.fits'):
            os.remove(filename)

def drizzle_f814w(clean_files=True):
    """ 

    """
    from drizzlepac import astrodrizzle
    import glob



    main_dir = '/home/egentry/Data/HST/PictorA/PSFs/f814w/'  

    artificial_star_dir = main_dir + 'artificial_stars/'
    drizzled_dir        = main_dir + 'drizzled/'
    psf_rootname        = 'psf'

    if not os.path.exists(drizzled_dir):
        os.makedirs(drizzled_dir)

    

    astrodrizzle.AstroDrizzle(artificial_star_dir + '*flt.fits', output=drizzled_dir + psf_rootname, runfile=drizzled_dir + 'astrodrizzle.log',
        preserve=False, static=False, skysub=True, driz_separate=False, median=False, blot=False, driz_cr=False, 
        driz_combine=True, final_wht_type='EXP', final_pixfrac=0.7, final_bits=4096,
        final_wcs=True, final_scale=0.01981)

    if clean_files == True:
        for filename in glob.glob(artificial_star_dir +'*_final_mask.fits'):
            os.remove(filename)
                        
                        


def crop_f160w():
    """
        Takes a drizzled TinyTim PSF, filters out NaN values, crops to 30arcsec box suitable for GALFIT
    """
    
    import numpy as np

    main_dir = '/home/egentry/Data/HST/PictorA/PSFs/drizzler/'

    artificial_star_dir = main_dir + 'artificial_stars/'
    drizzled_dir        = main_dir + 'drizzled/'
    cropped_dir         = main_dir + 'cropped/'
    psf_rootname        = 'psf'

    if not os.path.exists(cropped_dir):
        os.makedirs(cropped_dir)

    main_dir = '/home/egentry/Data/HST/PictorA/PSFs/f475w/'  

    flt_dir  = '/home/egentry/Data/HST/PictorA/IBJX01020/'

    artificial_star_dir = main_dir + 'artificial_stars/'
    drizzled_dir        = main_dir + 'drizzled/'
    psf_rootname        = 'psf'



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

    print valid_coords

    psf_data_cropped = psf_data[min_nonzero_y : max_nonzero_y + 1, min_nonzero_x : max_nonzero_x + 1]

    psf_cropped_hdulist = fits.PrimaryHDU(psf_data_cropped)
    psf_cropped_hdulist.writeto(cropped_dir + psf_filename + '_cropped.fits', clobber=True)
    
    psf_hdulist.close()


def crop_f475w():
    """
        Takes a drizzled TinyTim PSF, filters out NaN values, crops to 30arcsec box suitable for GALFIT
    """
    
    import numpy as np

    main_dir = '/home/egentry/Data/HST/PictorA/PSFs/f475w/'

    artificial_star_dir = main_dir + 'artificial_stars/'
    drizzled_dir        = main_dir + 'drizzled/'
    cropped_dir         = main_dir + 'cropped/'
    psf_rootname        = 'psf'

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

    valid_coords =  np.nonzero(psf_data[1:4000][1:4000])


    min_nonzero_x = np.amin( valid_coords[1] )
    max_nonzero_x = np.amax( valid_coords[1] )
    min_nonzero_y = np.amin( valid_coords[0] )
    max_nonzero_y = np.amax( valid_coords[0] )

    print valid_coords

    psf_data_cropped = psf_data[min_nonzero_y : max_nonzero_y + 1, min_nonzero_x : max_nonzero_x + 1]

    psf_cropped_hdulist = fits.PrimaryHDU(psf_data_cropped)
    psf_cropped_hdulist.writeto(cropped_dir + 'f475w_' + psf_filename + '_cropped.fits', clobber=True)
    
    psf_hdulist.close()

def crop_f814w():
    """
        Takes a drizzled TinyTim PSF, filters out NaN values, crops to 30arcsec box suitable for GALFIT
    """
    
    import numpy as np

    main_dir = '/home/egentry/Data/HST/PictorA/PSFs/f814w/'  

    artificial_star_dir = main_dir + 'artificial_stars/'
    drizzled_dir        = main_dir + 'drizzled/'
    cropped_dir         = main_dir + 'cropped/'
    psf_rootname        = 'psf'

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

    valid_coords =  np.nonzero(psf_data[1:4000][1:4000])


    min_nonzero_x = np.amin( valid_coords[1] )
    max_nonzero_x = np.amax( valid_coords[1] )
    min_nonzero_y = np.amin( valid_coords[0] )
    max_nonzero_y = np.amax( valid_coords[0] )

    print valid_coords

    psf_data_cropped = psf_data[min_nonzero_y : max_nonzero_y + 1, min_nonzero_x : max_nonzero_x + 1]

    psf_cropped_hdulist = fits.PrimaryHDU(psf_data_cropped)
    psf_cropped_hdulist.writeto(cropped_dir + 'f814w_' + psf_filename + '_cropped.fits', clobber=True)
    
    psf_hdulist.close()





def main_f160w():

    prepare_f160w()
    drizzle_f160w()
    crop_f160w()

def main_f475w():
    prepare_f475w()
    drizzle_f475w()
    crop_f475w()

def main_f814w():
    prepare_f814w()
    drizzle_f814w()
    crop_f814w()
    

def create_f475w():

    psf_dir = '/home/egentry/Data/HST/PictorA/PSFs/f475w/'

    agn_x = [1699, 1701, 1704]
    agn_y = [838, 841, 842]

    os.chdir(psf_dir)

    for (x,y) in izip(agn_x, agn_y):

        create_f475w_bash(x, y)

        subprocess.call('bash create_psf_f475w >/dev/null', shell=True) 
        subprocess.call('$TINYTIM/tiny2 params && $TINYTIM/tiny3 params', shell=True) 

        #clean up files:
        subprocess.call('rm *.tt3' ,shell=True)
        subprocess.call('rm *psf.fits' ,shell=True)

        subprocess.call('mv ' + str(x) + '_' + str(y) + '_00.fits' + ' ' + str(x) + '_' + str(y) + '.fits', shell=True) 

def create_f475w_bash(x,y):
    
    f = open('create_psf_f475w', 'w')

    f.write('#!/bin/bash' + '\n')
    f.write('' + '\n')
    f.write('$TINYTIM/tiny1 params <<EOF' + '\n')
    f.write('22' + '\n')
    f.write('2' + '\n')
    f.write(str(x) + ' ' + str(y) +  '\n')
    f.write('f475w' + '\n')
    f.write('3' + '\n')
    f.write('-1' + '\n')
    f.write('10' + '\n')
    f.write('0' + '\n')
    f.write(str(x) + '_' + str(y) + '_' + '\n') 
    f.write('EOF' + '\n')

    f.close()



def create_f814w():

    psf_dir = '/home/egentry/Data/HST/PictorA/PSFs/f814w/'

    agn_x = [1699, 1701, 1704]
    agn_y = [838, 841, 843]

    os.chdir(psf_dir)

    for (x,y) in izip(agn_x, agn_y):

        create_f814w_bash(x, y)

        subprocess.call('bash create_psf_f814w >/dev/null', shell=True) 
        subprocess.call('$TINYTIM/tiny2 params && $TINYTIM/tiny3 params', shell=True) 

        #clean up files:
        subprocess.call('rm *.tt3' ,shell=True)
        subprocess.call('rm *psf.fits' ,shell=True)

        subprocess.call('mv ' + str(x) + '_' + str(y) + '_00.fits' + ' ' + str(x) + '_' + str(y) + '.fits', shell=True) 

def create_f814w_bash(x,y):
    
    f = open('create_psf_f814w', 'w')

    f.write('#!/bin/bash' + '\n')
    f.write('' + '\n')
    f.write('$TINYTIM/tiny1 params <<EOF' + '\n')
    f.write('22' + '\n')
    f.write('2' + '\n')
    f.write(str(x) + ' ' + str(y) +  '\n')
    f.write('f475w' + '\n')
    f.write('3' + '\n')
    f.write('-1' + '\n')
    f.write('10' + '\n')
    f.write('0' + '\n')
    f.write(str(x) + '_' + str(y) + '_' + '\n') 
    f.write('EOF' + '\n')

    f.close()
    

#create_f475w()
#create_f814w()

main_f814w()
main_f475w()
main_f160w()

    
