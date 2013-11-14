main_dir = '/home/egentry/Data/HST/PictorA/IBJX01020/'
drizzled_dir = main_dir + 'drizzled/'


from astropy.io import fits
import os


def drizzle(clean_files=True):
    """ Drizzles FLT science images
        
        Currently using EXP weighting

        Version History
            2013.10.25 -- Eric Gentry -- Created script for F160w
            2013.11.14 -- Eric Gentry -- Adapted for F475W dataset

        To do:
            Combine dataset drizzlers into one script

    """

    from drizzlepac import astrodrizzle
    import glob

    if not os.path.exists(drizzled_dir):
        os.makedirs(drizzled_dir)

    os.chdir(main_dir) # required, because astrodrizzle can't handle path names with capital letters

    # 1st run
    astrodrizzle.AstroDrizzle('*flt.fits', output= drizzled_dir + 'f475w', driz_combine=False, clean=True)

    astrodrizzle.AstroDrizzle('*flt.fits', output= drizzled_dir + 'f475w', clean=True, driz_separate=False,  
        median=False, blot=False, driz_cr=False, driz_combine=True, final_wht_type='EXP',  
        final_pixfrac = 0.7, final_wcs=True, final_scale =0.06666 )

### Handled by using "clean=True" in astrodrizzle call:
#    if clean_files == True:
#        for filename in glob.glob(main_dir +'*_final_mask.fits'):
#            os.remove(filename)
#        for filename in glob.glob(main_dir +'tmp*.fits'):
#            os.remove(filename)
#        for filename in glob.glob(main_dir +'*crmask.fits'):
#            os.remove(filename)
#        for filename in glob.glob(main_dir +'*_blt.fits'):
#            os.remove(filename)
                        

def main():


    drizzle()

    #Filter error images?
    #Filter mask for inner/outer radii?

main()
   

    
