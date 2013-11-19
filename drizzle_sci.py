
from astropy.io import fits
from drizzlepac import astrodrizzle
import os


def drizzle(main_dir, drizzled_dir, filter_name, scale, clean_files=True):
    """ Drizzles FLT science images
        
        Currently using EXP weighting

        Version History
            2013.10.25 -- Eric Gentry -- Created script for F160w
            2013.11.14 -- Eric Gentry -- Adapted for F475W dataset
            2013.11.17 -- Eric Gentry -- Adapted for F814W dataset

        To do:
            Combine dataset drizzlers into one script

    """



    if not os.path.exists(drizzled_dir):
        os.makedirs(drizzled_dir)

    os.chdir(main_dir) # required, because astrodrizzle can't handle path names with capital letters

    # 1st run
    astrodrizzle.AstroDrizzle('*flt.fits', output= drizzled_dir + filter_name, driz_combine=False, clean=True, preserve=False)

    astrodrizzle.AstroDrizzle('*flt.fits', output= drizzled_dir + filter_name, clean=True, driz_separate=False, preserve=False, 
        median=False, blot=False, driz_cr=False, driz_combine=True, final_wht_type='EXP',  
        final_pixfrac = 0.7, final_wcs=True, final_scale = scale)
                        

def main_f160w():

    main_dir = '/home/egentry/Data/HST/PictorA/IBJX01010/'
    drizzled_dir = main_dir + 'drizzled/'
    filter_name = 'f160w'
    scale = 0.06666


    drizzle(main_dir, drizzled_dir, filter_name, scale)

def main_f475w():

    main_dir = '/home/egentry/Data/HST/PictorA/IBJX01020/'
    drizzled_dir = main_dir + 'drizzled/'
    filter_name = 'f475w'
    scale = .03962 / 2


    drizzle(main_dir, drizzled_dir, filter_name, scale)

def main_f814w():

    main_dir = '/home/egentry/Data/HST/PictorA/IBJX01030/'
    drizzled_dir = main_dir + 'drizzled/'
    filter_name = 'f814w'
    scale = .03962 / 2


    drizzle(main_dir, drizzled_dir, filter_name, scale)



def main():
    main_f160w()
    main_f475w()
    main_f814w()

main()
   

    
