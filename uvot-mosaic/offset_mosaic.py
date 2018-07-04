import numpy as np

from astropy.io import fits
from astropy import wcs
from astropy.stats import biweight_location, sigma_clip
from regions import read_ds9

import glob
import os
import os.path
import subprocess
import copy

from config_uvot_mosaic import __ROOT__


import pdb


def offset_mosaic(input_prefix,
                      output_prefix,
                      filter_list=['w2','m2','w1','uu','bb','vv'],
                      min_exp_w2=170, min_exp_m2=230, min_exp_w1=200,
                      min_exp_uu=0, min_exp_bb=0, min_exp_vv=0,
                      restack_id=False, mask_file=None):
    """
    Create mosaics in which the background varies between snapshots, so they need to be adjusted to match.

    The default minimum exposure times for the UV filters were estimated by LMZH to be reasonable values for M31, but no experimenting has been done with the optical filters, so their default is currently 0.

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

    restack_id : boolean (default = False)
        By default, if the stacking for a given target ID is done, it won't do the stacking again.  Set this to True to force the stacking to be done anyway.

    mask_file : string
        Name of ds9 region file with circles.  These areas will be masked when calculating the biweight of overlapping areas.

    Returns
    -------
    nothing

    """

    # make dictionary with the minimum exposure times
    min_exp = {'w2':min_exp_w2, 'm2':min_exp_m2, 'w1':min_exp_w1,
                   'uu':min_exp_uu, 'bb':min_exp_bb, 'vv':min_exp_vv}
    

    # go through each filter to build images

    for filt in filter_list:

        # ------------------------
        # find unique target IDs, and stack those first
        # ------------------------

        # open the images
        with fits.open(input_prefix + filt + '_sk_all.fits') as hdu_sk, fits.open(input_prefix + filt + '_ex_all.fits') as hdu_ex:

            # delete the 0th extensions (no images there, and they break later steps)
            del hdu_sk[0]
            del hdu_ex[0]
            
            # remove extensions with exposures shorter than minimum
            exp_time = np.array( [hdu_sk[i].header['EXPOSURE'] for i in range(len(hdu_sk))] )
            remove_ind = np.where(exp_time < min_exp[filt])[0]
            for ind in sorted(remove_ind, reverse=True):
                del hdu_sk[ind]
                del hdu_ex[ind]



            # all of the target IDs
            target_ids = np.array( [hdu_sk[i].header['TARG_ID'] for i in range(len(hdu_sk))] )
            # chop it down to just the unique ones
            target_ids = np.unique(target_ids)

           
            for targ in target_ids:

                print('')
                print('##### stacking target ID ' + str(targ) + ', filter ' + filt + ' #####')
                print('')

                # prefix for saving the files for this target ID
                file_prefix = output_prefix + str(targ) + '_' + filt

                # check if this one is done already (by looking for a count rate image)
                if os.path.isfile(file_prefix + '_cr.fits') and (restack_id == False):
                    print(str(targ)+' is already done')
                    print('')
                    continue
                
                
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
                if mask_file is not None:
                    temp_hdu_ex_adj = mask_image(temp_hdu_ex_adj, mask_file)

                # write out to files
                temp_hdu_sk.writeto('targ_temp_sk.fits', overwrite=True)
                temp_hdu_ex_adj.writeto('targ_temp_ex.fits', overwrite=True)
                
                # find the coordinates of the overlapping area
                overlap_x, overlap_y = find_overlap('targ_temp_ex.fits')

                # find the biweight of the overlapping areas
                biweight_cps = calc_overlap_val(temp_hdu_sk, temp_hdu_ex, overlap_x, overlap_y)

                # apply to the counts images
                hdu_sk_corr, _ = correct_sk(temp_hdu_sk, temp_hdu_ex, biweight_cps)

                # write out to files
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


        # output file names
        output_file_sk = output_prefix + filt + '_sk.fits'
        output_file_sk_all = output_prefix + filt + '_sk_all.fits'
        output_file_ex = output_prefix + filt + '_ex.fits'
        output_file_ex_all = output_prefix + filt + '_ex_all.fits'
        output_file_cr = output_prefix + filt + '_cr.fits'

        # start out the stacking with the first target ID
        subprocess.run('cp '+ output_prefix + str(target_ids[0]) + '_' + filt + '_sk.fits ' + output_file_sk, shell=True)
        subprocess.run('cp '+ output_prefix + str(target_ids[0]) + '_' + filt + '_ex.fits ' + output_file_ex, shell=True)
        subprocess.run('cp '+ output_prefix + str(target_ids[0]) + '_' + filt + '_sk.fits ' + output_file_sk_all, shell=True)
        subprocess.run('cp '+ output_prefix + str(target_ids[0]) + '_' + filt + '_ex.fits ' + output_file_ex_all, shell=True)
        # make a count rate image too
        with fits.open(output_file_sk) as h_sk, fits.open(output_file_ex) as h_ex:
            cr_hdu = fits.PrimaryHDU(data=h_sk[1].data/h_ex[1].data, header=h_sk[1].header)
            cr_hdu.writeto(output_file_cr, overwrite=True)

            
        # keep track of which target IDs still need to be appended to the image
        remaining_ids = copy.copy(target_ids[1:])


        # keep going while there are still IDs to append
        while len(remaining_ids) > 0:

            # file names for the target IDs
            remaining_id_files_sk = [output_prefix + str(t) + '_' + filt + '_sk.fits' for t in remaining_ids]
            remaining_id_files_ex = [output_prefix + str(t) + '_' + filt + '_ex.fits' for t in remaining_ids]
            
            # find the target ID that has the best overlap with current mosaic
            # (returns index and the overlapping pixels)
            best_ind, overlap_x, overlap_y = most_overlap(output_file_ex, remaining_id_files_ex)

            # make an HDU with the counts (sk) image for the mosaic and best ID
            with fits.open(output_file_sk) as hdu_mosaic_sk, fits.open(remaining_id_files_sk[best_ind]) as hdu_best_sk:
                temp_hdu_sk = fits.HDUList()
                temp_hdu_sk.append(fits.ImageHDU(data=hdu_mosaic_sk[1].data, header=hdu_mosaic_sk[1].header))
                temp_hdu_sk.append(fits.ImageHDU(data=hdu_best_sk[1].data, header=hdu_best_sk[1].header))
            # make an HDU with the exposure image for the mosaic and best ID
            with fits.open(output_file_ex) as hdu_mosaic_ex, fits.open(remaining_id_files_ex[best_ind]) as hdu_best_ex:
                temp_hdu_ex = fits.HDUList()
                temp_hdu_ex.append(fits.ImageHDU(data=hdu_mosaic_ex[1].data, header=hdu_mosaic_ex[1].header))
                temp_hdu_ex.append(fits.ImageHDU(data=hdu_best_ex[1].data, header=hdu_best_ex[1].header))
                 
            # find the biweight of the overlapping areas
            biweight_cps = calc_overlap_val(temp_hdu_sk, temp_hdu_ex, overlap_x, overlap_y)

            # apply to the counts images
            hdu_sk_corr, delta_cps = correct_sk(temp_hdu_sk, temp_hdu_ex, biweight_cps)

            # save those changes to the individual target ID segments
            with fits.open(output_file_sk_all) as hdu_sk_all, fits.open(output_file_ex_all) as hdu_ex_all:
                # adjust segments
                for h in range(1,len(hdu_sk_all)):
                    hdu_sk_all[h].data = (hdu_sk_all[h].data/hdu_ex_all[h].data + delta_cps[0]) * hdu_ex_all[h].data
                    hdu_sk_all[h].data[hdu_ex_all[h].data == 0] = 0
                hdu_sk_all.append(fits.ImageHDU(data=hdu_sk_corr[1].data, header=hdu_sk_corr[1].header))
                hdu_ex_all.append(fits.ImageHDU(data=temp_hdu_ex[1].data, header=temp_hdu_ex[1].header))
                # write out to files
                hdu_sk_all.writeto(output_file_sk_all, overwrite=True)
                hdu_ex_all.writeto(output_file_ex_all, overwrite=True)
                
                           
            # stack with uvotimsum
            cmd = 'uvotimsum ' + output_file_sk_all + ' ' + output_file_sk + ' exclude=none clobber=yes'
            subprocess.run(cmd, shell=True)
            cmd = 'uvotimsum ' + output_file_ex_all + ' ' + output_file_ex + ' method=EXPMAP exclude=none clobber=yes'
            subprocess.run(cmd, shell=True)

            # make a count rate image too
            with fits.open(output_file_sk) as h_sk, fits.open(output_file_ex) as h_ex:
                cr_hdu = fits.PrimaryHDU(data=h_sk[1].data/h_ex[1].data, header=h_sk[1].header)
                cr_hdu.writeto(output_file_cr, overwrite=True)

               
            # finally, remove this index from the remaining IDs list
            remaining_ids = np.delete(remaining_ids, best_ind)

            
      

def most_overlap(mosaic_ex, id_list, mask_file=None):
    """
    Find the index of the target ID that has the best overlap with current mosaic

    Parameters
    ----------
    mosaic_ex : string
        name of the fits file with the exposure map for the current mosaic

    id_list : string
        exposure map filenames for the target IDs that need to be checked against the mosaic

    mask_file : string
        Name of ds9 region file with circles.  These areas will be masked when calculating the biweight of overlapping areas.


    Returns
    -------
    best_ind : integer
        index of the target ID with the most overlap

    overlap_x, overlap_y : lists of float arrays
        for the best overlap, this is the output from find_overlap

    """

    # initialize arrays
    best_ind = -99
    overlap_x = [[],[]]
    overlap_y = [[],[]]
    

    # read in the mosaic exposure map
    with fits.open(mosaic_ex) as hdu_mosaic:

        # go through the target IDs
        for t in range(len(id_list)):

            # open the current exposure map
            with fits.open(id_list[t]) as hdu_id:

                # copy over the two extensions
                temp_hdu_ex = fits.HDUList()
                temp_hdu_ex.append(fits.ImageHDU(data=hdu_mosaic[1].data, header=hdu_mosaic[1].header))
                temp_hdu_ex.append(fits.ImageHDU(data=hdu_id[1].data, header=hdu_id[1].header))

                # turn it into 0s/1s
                temp_hdu_ex = exp_to_ones(temp_hdu_ex)

                # apply mask
                if mask_file is not None:
                    temp_hdu_ex = mask_image(temp_hdu_ex, mask_file)

                # write it out
                temp_hdu_ex.writeto('temp_overlap_ex.fits', overwrite=True)
                
                # find the overlap
                current_overlap_x, current_overlap_y = find_overlap('temp_overlap_ex.fits')

                # if there's no overlap, go to the next image
                if len(current_overlap_x) == 0:
                    continue

                # if there is overlap, and it's more than the biggest so far, save it
                if len(current_overlap_x[0]) > len(overlap_x[0]):
                    overlap_x = copy.deepcopy(current_overlap_x)
                    overlap_y = copy.deepcopy(current_overlap_y)
                    best_ind = copy.copy(t)
    
    # delete temp files
    subprocess.run('rm temp_overlap_ex.fits', shell=True)

    
    return best_ind, overlap_x, overlap_y



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
    Using a multi-extension fits file, in which each extension is 0 in places where there is no exposure and 1 where there is exposure, find a common overlapping area between ALL of the extensions.  If there is no area with 100% overlap, this will return empty lists.

    Parameters
    ----------
    image_footprint : string
        Name of the multi-extension fits file for which the data array is 1 where there's image data, and 0 otherwise

    Returns
    -------
    overlap_x, overlap_y : lists of float arrays
        The X/Y coordinates of the pixels where there is full overlap.  Each element of the list contains the array of values for each extension.  If there is no overlap, both have length of 0.
   
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
        y_stack, x_stack = np.where(hdu_stack[1].data > 0.9999*len(hdu_ex))

        # if there's no overlap, return empty lists
        if len(x_stack) == 0:
            return [], []

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
    



def calc_overlap_val(hdu_sk, hdu_ex, overlap_x, overlap_y, method='biweight'):
    """
    Find a representative value (for now, assuming biweight) for the overlapping area 

    Parameters
    ----------
    hdu_sk : astropy hdu object
        An HDU with the sky (counts) images

    hdu_ex : astropy hdu object
        An HDU with the exposure images

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
    
    for h in range(len(hdu_sk)):

        count_rate = hdu_sk[h].data / hdu_ex[h].data
        
        grab_pix = []

        # get the pixel values
        for i in range(len(overlap_x[h])):
            grab_pix.append(count_rate[ int(overlap_y[h][i]), int(overlap_x[h][i]) ])

        # do a sigma clip
        pix_clip = sigma_clip(np.array(grab_pix), sigma=2.5, iters=3)
            
        # calculate biweight
        #biweight_noclip = biweight_location(np.array(grab_pix))
        biweight_clip = biweight_location(pix_clip.data[~pix_clip.mask])
        
        val.append(biweight_clip)


    
    return np.array(val)

    

def correct_sk(hdu_sk, hdu_ex, overlap_cps):
    """
    We want all overlapping areas to have the same counts/sec, so adjust the counts image (using the known exposure time) accordingly

    Parameters
    ----------
    hdu_sk : astropy hdu object
        An HDU with the sky (counts) images

    hdu_ex : astropy hdu object
        An HDU with the exposure images

    overlap_cps : array of floats
        The biweight (or other calculation) counts/sec/pixel in the overlapping region

    Returns
    -------
    hdu_corr : astropy hdu object
        the same hdu as the input sky image, but with an offset applied

    delta_cps : array of floats
        the amount (in counts/sec/pixel) that each image was shifted

    """

    print('')
    print('  ** adjusting count images')
    print('')


    # we want everything to match the minimum counts/sec
    min_cps = np.min(overlap_cps)

    # array to keep track of how much each counts image is adjusted
    delta_cps = np.zeros(len(hdu_sk))
    
    
    for h in range(len(hdu_sk)):
        # do the offset
        delta_cps[h] = - (overlap_cps[h] - min_cps)
        hdu_sk[h].data = (hdu_sk[h].data/hdu_ex[h].data + delta_cps[h]) * hdu_ex[h].data
        # use exposure map to set border to 0
        hdu_sk[h].data[hdu_ex[h].data == 0] = 0
        
        

    return hdu_sk, delta_cps


        

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
