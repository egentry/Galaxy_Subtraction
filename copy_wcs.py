from astropy.io import fits
from astropy    import wcs
import numpy as np
import os


def copy_galfit_output_wcs(target_filename):
    """ Copies the WCS information to all frames of a galfit output file

        Version History
            2013.11.05 -- Eric Gentry -- Created script -- simply copied entire header over output image
            2013.11.17 -- Eric Gentry -- Script only copies specified FITS header cards

        Things to do:
            Create a function for copying wcs between arbitrary files with arbitrary output filename
            Check what WCS names are actually necessary / what we can or can't expect

    """
    #cards_to_copy was copied by hand from f160w Pictor A data that seemed relevant
    card_names_to_copy =    [   'WCSAXES',       
                                'CRPIX1',   
                                'CRPIX2',  
                                'CRVAL1',          
                                'CRVAL2',  
                                'CTYPE1',          
                                'CTYPE2', 
                                'ORIENTAT',
                                'VAFACTOR',
                                'CD1_1', 
                                'CD1_2',
                                'CD2_1',   
                                'CD2_2',   
                                'LTV1',        
                                'LTV2',            
                                'LTM1_1',     
                                'LTM2_2',   
                                'PA_APER', 
                                'RA_APER',
                                'DEC_APER',
                                'NCOMBINE',
                                'CENTERA1',
                                'CENTERA2',
                                'SIZAXIS1',
                                'SIZAXIS2',
                                'BINAXIS1',
                                'BINAXIS2',
                                'SAMPNUM',
                                'SAMPTIME',
                                'DELTATIM',
                                'ROUTTIME',
                                'TDFTRANS',
                                'PODPSFF',
                                'STDCFFF', 
                                'STDCFFP',
                                'WCSNAMEO',
                                'WCSAXESO',
                                'CRPIX1O',
                                'CRPIX2O',
                                'CDELT1O',
                                'CDELT2O', 
                                'CUNIT1O', 
                                'CUNIT2O', 
                                'CTYPE1O',
                                'CTYPE2O',
                                'CRVAL1O', 
                                'CRVAL2O',
                                'LONPOLEO',
                                'LATPOLEO',
                                'RESTFRQO',
                                'RESTWAVO',
                                'CD1_1O',
                                'CD1_2O', 
                                'CD2_1O', 
                                'CD2_2O', 
                                'WCSNAME'  ]


    overall_hdulist = fits.open(main_dir + target_filename)

    original_header = overall_hdulist[1].header
    model_header    = overall_hdulist[2].header
    resid_header    = overall_hdulist[3].header

    card_list = original_header.ascardlist()

    for card_name in card_names_to_copy:
        if card_name not in original_header:
            print 'Card not found in header: ', card_name
            print 'Check if it is needed?'
            return
        else:
#            print 'Found card: ', card_name
            card = card_list[card_name]
            model_header[card.key] = (card.value, card.comment)
            resid_header[card.key] = (card.value, card.comment)
            


    overall_hdulist.writeto(main_dir + target_filename, output_verify='fix', clobber=True)

    overall_hdulist.close()


def test():

    main_dir = '/home/egentry/Data/HST/PictorA/main_galaxy/tmp/full_sized/'

    filename_list = [   'tmp_galaxy_3.fits',
                        'tmp_galaxy_2.fits'    ]

#    for filename in filename_list:
#        copy(filename)

    copy_galfit_output_wcs(main_dir + filename_list[0])

test()

    
