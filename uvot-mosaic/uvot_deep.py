import numpy as np
from astropy.io import fits
import glob
import os
import subprocess

from config_uvot_mosaic import __ROOT__


import pdb



def uvot_deep(input_folders,
                  output_prefix,
                  filter_list=['w2','m2','w1','uu','bb','vv']):
    """
    For a set of UVOT images downloaded from HEASARC, do processing on each snapshot:
    * create a counts image
    * create exposure map
    * create LSS image
    * create mask image for bad pixels
    * create scattered light image

    Cases where an image will be skipped:
    * If imaging for a filter doesn't exist, it will be skipped, even if that filter name is in the input.
    * UVOT images are generally 2x2 binned.  If any images are unbinned, they will be skipped.
    * If a particular snapshot has no aspect correction, the astrometry is unreliable, so it will be skipped.
    

    Modeled off of Michael Siegel's code uvot_deep.pro


    Parameters
    ----------
    input_folders : list of strings
        each item of the string is the 11-digit name of the folder downloaded from HEASARC

    output_prefix : string
        the prefix for output files (be sure to include an underscore or similar for readability)

    filter_list : list of strings
        some or all of ['w2','m2','w1','uu','bb','vv'] (default is all of them)

    Returns
    -------
    nothing

    """

    # full path to most recent teldef files
    caldb = os.environ['CALDB']
    teldef = {'uu':sorted(glob.glob(caldb+'/data/swift/uvota/bcf/teldef/*uu*'))[-1],
                  'bb':sorted(glob.glob(caldb+'/data/swift/uvota/bcf/teldef/*bb*'))[-1],
                  'vv':sorted(glob.glob(caldb+'/data/swift/uvota/bcf/teldef/*vv*'))[-1],
                  'w1':sorted(glob.glob(caldb+'/data/swift/uvota/bcf/teldef/*w1*'))[-1],
                  'm2':sorted(glob.glob(caldb+'/data/swift/uvota/bcf/teldef/*m2*'))[-1],
                  'w2':sorted(glob.glob(caldb+'/data/swift/uvota/bcf/teldef/*w2*'))[-1] }
                  

    # ------------------------
    # identify the filters in each snapshot
    # ------------------------

    # dictionary to hold filters that exist for each folder
    filter_exist = {key:[] for key in input_folders}
    
    for i in input_folders:

        # list all of the sky images
        sk_list = glob.glob(i + '/uvot/image/*_sk.img')

        # grab the filter from the filename of each sky image
        for sk in sk_list:
            filter_name = sk[-9:-7]
            if filter_name in filter_list:
                filter_exist[i].append(filter_name)


    # ------------------------
    # go through each filter and build the images
    # ------------------------

    for filt in filter_list:

        # get the images that have observations in that filter
        obs_list = [im for im in filter_exist.keys() if filt in filter_exist[im]]

        for obs in obs_list:

            # counts image (labeled as sk)
            sk_image = obs+'/uvot/image/sw'+obs+'u'+filt+'_sk.img'
            # exposure image
            ex_image = obs+'/uvot/image/sw'+obs+'u'+filt+'_ex.img'
            # attitude file
            att_sat = obs+'/auxil/sw'+obs+'sat.fits'
            # make a more accurate attitude file
            att_uat = obs+'/auxil/sw'+obs+'uat.fits'
            corr_file = obs+'/uvot/hk/sw'+obs+'uac.hk'
            if not os.path.isfile(att_uat):
                cmd = 'uvotattcorr attfile=' + att_sat + ' corrfile=' + corr_file + ' outfile=' + att_uat
                subprocess.run(cmd, shell=True)


            # scattered light images
            scattered_light(obs, filt, teldef[filt])

            # mask images (including masking the exposure maps)
            mask_image(obs, filt, teldef[filt])
            
            # LSS images

            
            pdb.set_trace()



def scattered_light(obs_folder, obs_filter, teldef_file):

    """
    Create scattered light images with the same orientation as the input snapshots

    Parameters
    ----------
    obs_folder : string
        the 11-digit name of the folder downloaded from HEASARC

    obs_filter : string
        one of the UVOT filters ['w2','m2','w1','uu','bb','vv']
    
    teldef_file : string
        full path+name for the teldef file


    Returns
    -------
    nothing

    """

    # counts image (labeled as sk)
    sk_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_sk.img'
    # attitude files
    att_uat = obs_folder+'/auxil/sw'+obs_folder+'uat.fits'
    att_sat = obs_folder+'/auxil/sw'+obs_folder+'sat.fits'

    
    with fits.open(sk_image) as hdu_sk:
        
        # grab the ra/dec/roll from the sky file
        ra_pnt = str(hdu_sk[0].header['RA_PNT'])
        dec_pnt = str(hdu_sk[0].header['DEC_PNT'])
        roll_pnt = str(hdu_sk[0].header['PA_PNT'])

        # create HDU for the scattered light images
        hdu_sl = fits.HDUList()
        # copy over the primary header from the sky image
        hdu_sl.append(fits.PrimaryHDU(header=hdu_sk[0].header))

        # for each image extension, make the SL image
        for i in range(1,len(hdu_sk)):
 
            # create image
            skytime = '{:.7f}'.format( (hdu_sk[i].header['TSTART'] + hdu_sk[i].header['TSTOP'])/2 )
            cmd = 'swiftxform infile='+__ROOT__+'/scattered_light_images/scal_'+obs_filter+'_smooth_2x2.fits' + \
                  ' outfile=temp.sl attfile=' + att_sat + ' teldeffile=' + teldef_file + ' method=AREA' + \
                  ' to=sky clobber=yes bitpix=-32 ra='+ra_pnt + ' dec='+dec_pnt + ' roll='+roll_pnt + \
                  ' skytime=MET:'+skytime
            subprocess.run(cmd, shell=True)
            
            # append it to the big fits file
            with fits.open('temp.sl') as hdu_sl_slice:
                hdu_sl.append(fits.ImageHDU(data=hdu_sl_slice[0].data, header=hdu_sl_slice[0].header))
                
            # delete the image
            os.remove('temp.sl')
            

        # write out all of the compiled SL images
        hdu_sl.writeto(obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'.sl', overwrite=True)

    subprocess.run('rm .nfs*', shell=True)
    




def mask_image(obs_folder, obs_filter, teldef_file):

    """
    Create a bad pixel map, and use that to create mask images and masked exposure maps

    Parameters
    ----------
    obs_folder : string
        the 11-digit name of the folder downloaded from HEASARC

    obs_filter : string
        one of the UVOT filters ['w2','m2','w1','uu','bb','vv']
    
    teldef_file : string
        full path+name for the teldef file


    Returns
    -------
    nothing

    """

    # counts image (labeled as sk)
    sk_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_sk.img'
    #exposure image
    ex_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_ex.img'
    # attitude files
    att_uat = obs_folder+'/auxil/sw'+obs_folder+'uat.fits'
    att_sat = obs_folder+'/auxil/sw'+obs_folder+'sat.fits'


    # make bad pixel map
    bad_pix = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'.badpix'
    subprocess.run('uvotbadpix infile='+sk_image + ' badpixlist=CALDB' + \
                   ' outfile='+bad_pix + ' clobber=yes', shell=True)

    # regenerate exposure maps
    ex_image_new = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_ex_mask.img'
    mask_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_mask.img'
    cmd = 'uvotexpmap infile='+sk_image + ' outfile='+ex_image_new + ' maskfile='+mask_image + \
          ' badpixfile='+bad_pix + ' method=MEANFOV attfile='+att_sat + ' teldeffile='+teldef_file + \
          ' masktrim=25 clobber=yes'




    
def stuff(sk_image, ex_image, att_file, teldef_file):

    """
    Create a scattered light image with the same orientation as the input snapshot

    Parameters
    ----------
    sk_image : string
        path+name for the sky (counts) image

    ex_image : string
        path+name for the exposure map image

    att_file : string
        path+name for the attitude file

    teldef_file : string
        full path+name for the teldef file


    Returns
    -------
    image : array
        the scattered light image

    header : something
        the header associated with that image
    
    """

if __name__ == '__main__':

    uvot_deep(['00037723001','00037723002'], 'test_', ['w2','m2','w1'])
