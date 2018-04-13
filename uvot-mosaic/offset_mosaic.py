import numpy as np

from astropy.io import fits
from astropy import wcs
from astropy.stats import biweight_location, sigma_clip
from regions import read_ds9

import glob
import os
import subprocess
import copy

from config_uvot_mosaic import __ROOT__


import pdb


def offset_mosaic(input_prefix,
                      output_prefix,
                      filter_list=['w2','m2','w1','uu','bb','vv'],
                      min_exp_w2=170, min_exp_m2=230, min_exp_w1=200,
                      min_exp_uu=0, min_exp_bb=0, min_exp_vv=0):
    """
    Create mosaics in which the background varies between snapshots, so they need to be adjusted to match.

    The default minimum exposure times for the UV filters were estimated by LMZH to be reasonable values, but no experimenting has been done with the optical filters, so their default is currently 0.

    Parameters
    ----------
    input_prefix : string
        The prefix for the input images (the multi-extension fits files created by uvot_deep)

    output_prefix : string
        The prefix for output files (be sure to include an underscore or similar for readability)

    filter_list : list of strings
        Some or all of ['w2','m2','w1','uu','bb','vv'] (default is all of them)

    min_exp_w2, min_exp_m2, min_exp_w1, min_exp_uu, min_exp_bb, min_exp_vv : integers
        Minimum exposure times (in seconds) for each filter.  Any snapshots with exposure times shorter than this will be discarded.

    Returns
    -------
    nothing

    """

    # make dictionary with the minimum exposure times
    min_exp = {'w2':min_exp_w2, 'm2':min_exp_m2, 'w1':min_exp_w1,
                   'uu':min_exp_uu, 'bb':min_exp_bb, 'vv':min_exp_vv}
    

    # go through each filter to build images

    for filt in filter_list:

        # open the images
        with fits.open(input_prefix + filt + '_sk_all.fits') as hdu_sk, fits.open(input_prefix + filt + '_ex_all.fits') as hdu_ex:

            # remove extensions with exposures shorter than minimum
            exp_time = np.array( [hdu_sk[i].header['EXPOSURE'] for i in range(len(hdu_sk))] )
            remove_ind = np.where(exp_time < min_exp[filt])[0]
            for ind in sorted(remove_ind, reverse=True):
                del hdu_sk[ind]
                del hdu_ex[ind]


            # ------------------------
            # find unique target IDs, and stack those first
            # ------------------------

            # all of the target IDs
            target_id = np.array( [hdu_sk[i].header['TARG_ID'] for i in range(len(hdu_sk))] )

           
            for targ in np.unique(target_id):
                
                # temp file to hold snapshots with current target ID
                temp_hdu_sk = fits.HDUList()
                temp_hdu_ex = fits.HDUList()

                # append matching snapshots
                [temp_hdu_sk.append(fits.ImageHDU(data=hdu_sk[i].data, header=hdu_sk[i].header)) for i in range(len(hdu_sk)) if hdu_sk[i].header['TARG_ID'] == targ]
                [temp_hdu_ex.append(fits.ImageHDU(data=hdu_ex[i].data, header=hdu_ex[i].header)) for i in range(len(hdu_sk)) if hdu_sk[i].header['TARG_ID'] == targ]

                # turn exposure maps into 0s and 1s
                temp_hdu_ex_adj = copy.deepcopy(temp_hdu_ex)
                temp_hdu_ex_adj = exp_to_ones(temp_hdu_ex_adj)

                # mask areas with foreground stars, etc.
                temp_hdu_ex_adj = mask_image(temp_hdu_ex_adj, 'bg_mask.reg')

                # write out to files
                temp_hdu_sk.writeto('targ_temp_sk.fits', overwrite=True)
                temp_hdu_ex_adj.writeto('targ_temp_ex.fits', overwrite=True)
                
                # find the coordinates of the overlapping area
                overlap_x, overlap_y = find_overlap('targ_temp_ex.fits')

                # find the biweight of the overlapping areas
                biweight_counts = calc_overlap_val(temp_hdu_sk, overlap_x, overlap_y)
                print('biweights: ', biweight_counts)
                pdb.set_trace()

                # apply to the counts images
                hdu_sk_corr = correct_sk(temp_hdu_sk, biweight_counts, temp_hdu_ex)

                # write out to files
                file_prefix = output_prefix + str(targ) + '_' + filt
                hdu_sk_corr.writeto(file_prefix + '_sk_all.fits', overwrite=True)
                temp_hdu_ex.writeto(file_prefix + '_ex_all.fits', overwrite=True)
                
                # stack with uvotimsum
                cmd = 'uvotimsum ' + file_prefix + '_sk_all.fits ' + \
                    file_prefix + '_sk.fits exclude=none clobber=yes'
                subprocess.run(cmd, shell=True)
                cmd = 'uvotimsum ' + file_prefix + '_ex_all.fits ' + \
                    file_prefix + '_ex.fits method=EXPMAP exclude=none clobber=yes'
                subprocess.run(cmd, shell=True)

                # make a count rate image too
                with fits.open(file_prefix + '_sk.fits') as h_sk, fits.open(file_prefix + '_ex.fits') as h_ex:
                    cr_hdu = fits.PrimaryHDU(data=h_sk[1].data/h_ex[1].data, header=h_sk[1].header)
                    cr_hdu.writeto(file_prefix + '_cr.fits', overwrite=True)
               
                # delete temporary files
                subprocess.run('rm targ_temp_*.fits', shell=True)
                
               
            # ------------------------
            # combine the stacks
            # ------------------------

            

            
                
def mask_image(hdu, reg_file):
    """
    Mask images using a ds9 region file

    Parameters
    ----------
    hdu : astropy hdu object
        An HDU with however many extensions you want

    reg_file : string
        filename of a ds9 region file (all circles)


    Returns
    -------
    hdu_mask : astropy hdu object
        the same hdu as the input, but with 0s 

    """
    
    print('')
    print('  ** masking HDU extensions')
    print('')

    # read in the ds9 file
    regions = read_ds9(reg_file)
    # get ra/dec/radius (all in degrees)
    reg_ra = np.array( [regions[i].center.ra.deg[0] for i in range(len(regions))] )
    reg_dec = np.array( [regions[i].center.dec.deg[0] for i in range(len(regions))] )
    reg_rad_deg = np.array( [regions[i].radius.value for i in range(len(regions))] )/3600

    # go through each extension and mask
    for h in hdu:
        
        wcs_h = wcs.WCS(h.header)

        # convert to x/y
        reg_x, reg_y = wcs_h.wcs_world2pix(reg_ra, reg_dec, 1)
        reg_rad_pix = reg_rad_deg / wcs.utils.proj_plane_pixel_scales(wcs_h)[0]

        # to save time, do coarse removal of regions outside of image
        outside = (reg_x < -200) | (reg_x > h.data.shape[0]+200) | (reg_y < -200) | (reg_y > h.data.shape[1]+200)
        reg_x = reg_x[~outside]
        reg_y = reg_y[~outside]
        reg_rad_pix = reg_rad_pix[~outside]

        # go through the regions and mask
        height, width = h.data.shape
        y_grid, x_grid = np.ogrid[:height, :width]

        for i in range(len(reg_x)):
            dist_from_center = np.sqrt((x_grid - reg_x[i])**2 + (y_grid-reg_y[i])**2)
            mask = dist_from_center <= reg_rad_pix[i]
            h.data[mask] = 0

    return hdu



    
def find_overlap(image_footprint):
    """
    Using a multi-extension fits file, which is 0 in places where there is no exposure, find a common overlapping area.  This particular code is optimized for images that are assumed to have a common overlapping area.

    Parameters
    ----------
    image_footprint : string
        Name of the multi-extension fits file for which the data array is 1 where there's image data, and 0 otherwise

    Returns
    -------
    overlap_x, overlap_y : lists of float arrays
        The X/Y coordinates of the pixels where there is full overlap.  Each element of the list contains the array of values for each extension.
   
    """

    print('')
    print('  ** identifying overlapping area')
    print('')

    # add up the maps
    cmd = 'uvotimsum ' + image_footprint + ' targ_temp_ex_stack.fits method=EXPMAP exclude=none clobber=yes'
    subprocess.run(cmd, shell=True)

    # read in the stack and the original file
    with fits.open('targ_temp_ex_stack.fits') as hdu_stack, fits.open(image_footprint) as hdu_ex:

        # x/y coordinates of max overlap in stacked image
        y_stack, x_stack = np.where(hdu_stack[1].data > 0.999*len(hdu_ex))

        # convert to ra/dec
        wcs_stack = wcs.WCS(hdu_stack[1].header)
        ra_stack, dec_stack = wcs_stack.wcs_pix2world(x_stack, y_stack, 1)

        # for each footprint image, convert back to corresponding x/y
        overlap_x = []
        overlap_y = []
        
        for h in hdu_ex:
            wcs_ex = wcs.WCS(h.header)
            ov_x, ov_y = wcs_ex.wcs_world2pix(ra_stack, dec_stack, 1)

            overlap_x.append(ov_x)
            overlap_y.append(ov_y)

        
    # return the list of arrays
    return overlap_x, overlap_y
    



def calc_overlap_val(hdu, overlap_x, overlap_y, method='biweight'):
    """
    Find a representative value (for now, assuming biweight) for the overlapping area 

    Parameters
    ----------
    hdu : astropy hdu object
        An HDU with the sky (counts) images

    overlap_x, overlap_y : lists of float arrays
        The X/Y coordinates of the pixels where there is full overlap.  Each element of the list contains the array of values for each extension.

    method : string (default='biweight')
        The calculation to use - later I might add median or mean

    Returns
    -------
    val : array of floats
       The result of the calculation for each extension in the HDU

    """

    print('')
    print('  ** calculating overlap ' + method)
    print('')


    val = []
    
    for h in range(len(hdu)):

        grab_pix = []

        # get the pixel values
        for i in range(len(overlap_x[h])):
            grab_pix.append(hdu[h].data[ int(overlap_y[h][i]), int(overlap_x[h][i]) ])

        # do a sigma clip
        pix_clip = sigma_clip(np.array(grab_pix), sigma=2.5, iters=3)
            
        # calculate biweight
        #biweight_noclip = biweight_location(np.array(grab_pix))
        biweight_clip = biweight_location(pix_clip.data[~pix_clip.mask])
        
        val.append(biweight_clip)


    
    return np.array(val)

    

def correct_sk(hdu_sk, overlap_counts, hdu_ex):
    """
    We want all overlapping areas to have the same counts/sec, so adjust the counts image (using the known exposure time) accordingly

    Parameters
    ----------
    hdu_sk : astropy hdu object
        An HDU with the sky (counts) images

    overlap_counts : array of floats
        The average (or other calculation) counts/pixel in the overlapping region

    hdu_ex : astropy hdu object
        An HDU with the exposure images

    Returns
    -------
    hdu_corr : astropy hdu object
        the same hdu as the input, but with an offset applied

    """

    print('')
    print('  ** adjusting count images')
    print('')


    # get exposure times
    exp_time = np.array( [h.header['EXPOSURE'] for h in hdu_sk] )

    # convert to counts/sec
    overlap_cps = overlap_counts / exp_time

    # we want everything to match the minimum counts/sec
    min_cps = np.min(overlap_cps)
    
    
    for h in range(len(hdu_sk)):
        # do the offset
        hdu_sk[h].data = (hdu_sk[h].data/exp_time[h] - (overlap_cps[h] - min_cps)) * exp_time[h]
        # use exposure map to set border to 0
        hdu_sk[h].data[hdu_ex[h].data == 0] = 0
        
        

    return hdu_sk


        

def overlap_stack(hdu):
    """
    Starting with the first extension, find the extension that overlaps most, offset it, and stack
    """



def exp_to_ones(hdu):
    """
    Turn exposure maps into 1s rather than the exposure time

    Input is the HDU object (NOT the saved fits file)
    """

    for h in range(len(hdu)):
        exp = hdu[h].data
        exp[exp > 0] = 1

    return hdu
