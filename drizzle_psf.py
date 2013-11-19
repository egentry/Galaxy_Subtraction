
from astropy.io import fits
import os
import subprocess
from itertools import izip

base_dir = '/home/egentry/Data/HST/PictorA/PSFs/'

def prepare(flt_dir, filter_name, agn_x, agn_y, flt_indices):
    """ Prepare PSFs to be drizzled 

        Expects FLT + PSF fits files
            Number of FLT files should match number of PSF files
            FLT files should be multiextension .fits files, as given by NASA ('SCI' = index 1)

        Input pixels follow a one-indexing convention, to match DS9. Pixels used by numpy follow 0-indexing convention


        Version History
            2013.08.09 -- Eric Gentry -- Created limited purpose draft
            2013.11.19 -- Eric Gentry -- Generalized for all datasets of Pictor A

        To do before usign on other datasets:
            Change format of filename of flt files

    """

    main_dir = base_dir + filter_name + '/'  

    artificial_star_dir = main_dir + 'artificial_stars/'
    drizzled_dir        = main_dir + 'drizzled/'
    cropped_dir         = main_dir + 'cropped/'
    psf_rootname        = 'psf'

    if not os.path.exists(artificial_star_dir):
        os.makedirs(artificial_star_dir)


    for (x,y, i) in izip(agn_x, agn_y, flt_indices):
    
        psf_filename =         '{0}_{1}.fits'.format(x,y)
        flt_filename =  'ibjx01y{}q_flt.fits'.format(i)  # To be done later: accept more general filenames
        
        psf_hdulist = fits.open(main_dir + psf_filename)
        
        flt_hdulist = fits.open(flt_dir + flt_filename)

        
        flt_data = flt_hdulist[1].data
        psf_data = psf_hdulist[0].data

        flt_size_x = flt_data.shape[1]
        flt_size_y = flt_data.shape[0]

        psf_size_x = psf_data.shape[1]
        psf_size_y = psf_data.shape[0]

        psf_center_x = round(psf_size_x / 2.0)
        psf_center_y = round(psf_size_y / 2.0)

        if filter_name == 'f160w':
            flt_data.fill(0)    
        elif filter_name == 'f475w' or filter_name == 'f814w':
            flt_data.fill(0)
            flt_hdulist[3].data.fill(0)     #also fill second chip (CCDCHIP=1) with 0's
        else:
            print 'Unrecognized filter name: ', filter_name

        #0 indexed for numpy arrays:
        min_nonzero_x = x-1 - (psf_center_x -1)
        max_nonzero_x = x-1 + (psf_center_x -1) -1
        min_nonzero_y = y-1 - (psf_center_y -1)
        max_nonzero_y = y-1 + (psf_center_y -1) -1

        for x_i in xrange(psf_size_x):
            for y_i in xrange(psf_size_y):
                flt_data[y_i + min_nonzero_y][x_i + min_nonzero_x] = psf_data[y_i][x_i] * 1e6

        flt_hdulist.writeto(artificial_star_dir + flt_filename, clobber=True)
        flt_hdulist.close()
        psf_hdulist.close()



def drizzle(filter_name, scale, clean_files=True):
    """ 

    """

    from drizzlepac import astrodrizzle
    import glob



    main_dir = base_dir + filter_name + '/'  

    artificial_star_dir = main_dir + 'artificial_stars/'
    drizzled_dir        = main_dir + 'drizzled/'
    psf_rootname        = 'psf'

    if not os.path.exists(drizzled_dir):
        os.makedirs(drizzled_dir)

    

    astrodrizzle.AstroDrizzle(artificial_star_dir + '*flt.fits', output=drizzled_dir + 'psf', runfile=drizzled_dir + 'astrodrizzle.log',
        preserve=False, static=False, skysub=True, driz_separate=False, median=False, blot=False, driz_cr=False, 
        driz_combine=True, final_wht_type='EXP', final_pixfrac= 0.7, final_bits= 4096,
        final_wcs=True, final_scale=scale)

    if clean_files == True:
        for filename in glob.glob(artificial_star_dir +'*_final_mask.fits'):
            os.remove(filename)




def crop(filter_name, max_x, max_y):
    """
        Takes a drizzled TinyTim PSF, filters out NaN values, crops to 30arcsec box suitable for GALFIT
    """
    
    import numpy as np

    main_dir = base_dir + filter_name + '/'  

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

    valid_coords =  np.nonzero(psf_data[0:max_y, 0:max_x])



    min_nonzero_x = np.amin( valid_coords[1] )
    max_nonzero_x = np.amax( valid_coords[1] )
    min_nonzero_y = np.amin( valid_coords[0] )
    max_nonzero_y = np.amax( valid_coords[0] )

#    print valid_coords

    psf_data_cropped = psf_data[min_nonzero_y : max_nonzero_y + 1, min_nonzero_x : max_nonzero_x + 1]

    psf_cropped_hdulist = fits.PrimaryHDU(psf_data_cropped)
    psf_cropped_hdulist.writeto(cropped_dir + filter_name + '_' + psf_filename + '_cropped.fits', clobber=True)
    
    psf_hdulist.close()
 

def create_f160w():

    psf_dir = '/home/egentry/Data/HST/PictorA/PSFs/f160w/'

    agn_centers = [ (363, 143), (366, 146) , (369, 150)]

    agn_x = [364, 367, 370]
    agn_y = [144, 147, 151]

    os.chdir(psf_dir)

    for (x,y) in izip(agn_x, agn_y):

        create_f160w_bash(x, y)

        subprocess.call('bash create_psf_f160w >/dev/null', shell=True) 
        subprocess.call('$TINYTIM/tiny2 params && $TINYTIM/tiny3 params', shell=True) 

        #clean up files:
        subprocess.call('rm *.tt3' ,shell=True)
        subprocess.call('rm *psf.fits' ,shell=True)

        subprocess.call('mv ' + str(x) + '_' + str(y) + '_00.fits' + ' ' + str(x) + '_' + str(y) + '.fits', shell=True) 

def create_f160w_bash(x,y):
    
    f = open('create_psf_f160w', 'w')

    f.write('#!/bin/bash' + '\n')
    f.write('' + '\n')
    f.write('$TINYTIM/tiny1 params <<EOF' + '\n')
    f.write('23' + '\n')
    f.write(str(x) + ' ' + str(y) +  '\n')
    f.write('f160w' + '\n')
    f.write('3' + '\n')
    f.write('-1' + '\n')
    f.write('25' + '\n')
    f.write('0' + '\n')
    f.write(str(x) + '_' + str(y) + '_' + '\n') 
    f.write('EOF' + '\n')

    f.close()
   

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


def main_f160w():
    flt_dir  = '/home/egentry/Data/HST/PictorA/IBJX01010/'
    filter_name = 'f160w'

    agn_x = [364, 367, 370]
    agn_y = [144, 147, 151]

    flt_indices = [0, 1, 2]

    scale = 0.06666

    max_x = -1
    max_y = -1
    

    prepare(flt_dir, filter_name, agn_x, agn_y, flt_indices)
    drizzle(filter_name, scale)
    crop(filter_name, max_x, max_y)

def main_f475w():
    flt_dir  = '/home/egentry/Data/HST/PictorA/IBJX01020/'
    filter_name = 'f475w'

    agn_x = [1699, 1701, 1704]
    agn_y = [ 838,  841,  842]

    flt_indices = [4,6,8]

    scale = 0.01981

    max_x = 4000
    max_y = 4000
    

    prepare(flt_dir, filter_name, agn_x, agn_y, flt_indices)
    drizzle(filter_name, scale)
    crop(filter_name, max_x, max_y)


def main_f814w():
    flt_dir  = '/home/egentry/Data/HST/PictorA/IBJX01030/'
    filter_name = 'f814w'

    agn_x = [1699, 1701, 1704]
    agn_y = [ 838,  841,  843]

    flt_indices = ['a', 'c', 'e']

    scale = 0.01981

    max_x = 4000
    max_y = 4000
    

    prepare(flt_dir, filter_name, agn_x, agn_y, flt_indices)
    drizzle(filter_name, scale)
    crop(filter_name, max_x, max_y)
    
#create_f160w()
#create_f475w()
#create_f814w()

main_f814w()
main_f475w()
main_f160w()



    
