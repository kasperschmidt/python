# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
import numpy as np
import time
import pdb
import os, sys, inspect
from os import listdir
import copy
import shutil
from astropy.io import fits
from astropy import wcs
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import matplotlib.image as mpimg
from wrapper_lib import *
from wrapper_lib_psf import *
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

"""
Script to run the galfit wrapper GUI by Josie Kerutt.
Based on https://bitbucket.org/Leviosa/galfit_wrapper/src/master/wrapper_each_obj.py
but altered to read input from setup file instead of using hard-coded keywords.

--- INPUT ---
setupfile              Setup file specifying the individual inputs need for the GALFIT wrapper.

--- EXAMPLE OF USE ---

wrappersetup='/Path/to/setupfile/galfit_wrapper_setup_MyFavoriteObject.txt'
python /path/to/script/run_script_galfit_wrapper.py $wrappersetup

(or use alias:  galfit_wrapper_run $wrappersetup)

"""
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
from_1 = 0  ### CHANGE! If you want to start from specific object
to_1 = None # for everything write None
# - - - - - - - - - - - - - - - - - - -  HANDLE COMMAND LINE INPUTS - - - - - - - - - - - - - - - - - - -
Narg = len(sys.argv)-1
print(' - Found '+str(Narg)+' arguments:\n   '+str(sys.argv[1:]))
if Narg != 1:
    print('   ERROR: Expects just 1 argument, i.e., provide 1 wrapper setup file to run')
    sys.exit()
else:
    setupfile = sys.argv[1]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print('\n - Hence will run the stupfile:\n   '+str(setupfile))
setupinfo = np.genfromtxt(setupfile,dtype=None,names=True,comments='#',skip_header=1)
setupdic  = {}
for nn, name in enumerate(setupinfo['name']):
    setupdic[name] = setupinfo['parameter'][nn]

for key in setupdic.keys():
    if setupdic[key].lower() == 'none':
        setupdic[key] = None
    elif setupdic[key].lower() == 'true':
        setupdic[key] = True
    elif setupdic[key].lower() == 'false':
        setupdic[key] = False

use_for                 = 'notpsfrun' # setupdic['use_for']
field                   = setupdic['field']
band                    = setupdic['band']
area                    = setupdic['area']
input_catalogue         = setupdic['input_catalogue']
input_image             = setupdic['input_image']
sigma_image             = setupdic['sigma_image']
output_image            = setupdic['output_image']
Guo_catalogue           = setupdic['Guo_catalogue']
counterparts_catalogue  = setupdic['counterparts_catalogue']
counterparts_imgs       = setupdic['counterparts_imgs']
counterparts_suff       = setupdic['counterparts_suff']
output_catalogue        = setupdic['output_catalogue']
galfit_path             = setupdic['galfit_path']
image_dir               = setupdic['image_dir']
out_log_dir             = setupdic['out_log_dir']
result_folder           = setupdic['result_folder']
psf_image               = setupdic['psf_image']
psf_sampling            = setupdic['psf_sampling']
bad_pixel               = setupdic['bad_pixel']
param_constr            = setupdic['param_constr']
save_stuff              = setupdic['save_stuff']
results_final           = setupdic['results_final']
file_comments           = setupdic['file_comments']
use_shape_params        = setupdic['use_shape_params']
placeholder_bombed      = setupdic['placeholder_bombed']
input_params_dir        = setupdic['input_params_dir']
input_params            = setupdic['input_params']
try:
    input_phot_zero_ab  = setupdic['phot_zeropoint_ab']
except:
    input_phot_zero_ab  = None
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if 'margin' in setupdic.keys():
    margin      = float(setupdic['margin'])
    margin_unit = setupdic['margin_unit'].lower()
else:
    margin = 1.0 # to add to the area that is fit on both sides, in arcsec (change size of cutout)
    margin_unit = 'arcsec'
xconv  = 150. # Size of the convolution box (x y)
yconv  = 150.

# Sky
sky_value = -0.0000618255 # write NONE of you want no sky
sky_dx = 0.0
sky_dy = 0.0
fix_sky = False #False

disp_type = 'regular' # Display type (regular, curses, both)
options = 0 # Options: 0=normal run; 1,2=make model/imgblock & quit

fix_pos = False # write True or False
fit_others = False # write True of you want to fit other close objects

# write 'sersic', 'expdisk', 'gaussian', 'moffat' or 'nuker' or 'psf'
obj_type = 'sersic'
if use_for == 'psf':
    obj_type = 'psf'
    use_shape_params = None

# what to fit
ident = False #'Lya' # write False if you want to fit all objects

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if use_for == 'psf':
    psf_image = None
    obj_type = 'moffat'
    margin = 2.0

    print('   '+field+'_'+band)

    try:
        star_info = stars_field[field+'_'+band]
    except:
        raise ValueError('There is no PSF star for field '+field+'!')
        exit()
        
    ra_sn_psf = star_info['ra']
    dec_sn_psf = star_info['dec']
    psf_id = star_info['id']
    x_pos_psf = star_info['x']
    y_pos_psf = star_info['y']
    
    # You can also insert information for a new star here!

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### PART 1 - create an input parameter file

# read information from header

hdu_hst      = fits.open(input_image)
header_hst   = hdu_hst[0].header
data_hst     = hdu_hst[0].data

try:
    try:
        cdelt1   = header_hst['CDELT1']
        cdelt2   = header_hst['CDELT2']
    except:
        cdelt1   = header_hst['CD1_1']
        cdelt2   = header_hst['CD2_2']

    if cdelt1 == 1 or cdelt2 == 1:
        cdelt1   = header_hst['CD1_1']
        cdelt2   = header_hst['CD2_2']

    if input_phot_zero_ab is None:
        # Manually setting zeropoint if not in image FITS header
        if 'acs_775w_udf-03_cut_rot.fits' in input_image:
            phot_zero_ab = 25.69 # from https://archive.stsci.edu/prepds/hlf/#products (Hubble Legacy Fields)
        elif 'cutout_rings.v3.skycell' in input_image:
            phot_zero_ab = 25.00 # from "HIERARCH FPA.ZP" keyword in image header
        else:
            photzpt   = header_hst['PHOTZPT'] # ST magnitude photometric zeropoint
            photflam  = header_hst['PHOTFLAM']
            photplam  = header_hst['PHOTPLAM']
            phot_zero = -2.5*np.log10(photflam)+photzpt
            # convert st mag zeropoint to ab mag
            # this is taken from http://www.stsci.edu/hst/wfc3/phot_zp_lbn
            phot_zero_ab = -2.5*np.log10(photflam)-21.10-5*np.log10(photplam)+18.6921
    else:
        print('\n - Will use photometric zeropoint from setup file: phot_zeropoint_ab = '+str(input_phot_zero_ab))
        phot_zero_ab = input_phot_zero_ab
except:
    print('   There was no CDELT1, CDELT2 or PHOTZPT in the header of your file. ('+str(input_image)+')')
    import pdb; pdb.set_trace()
    exit()
hdu_hst.close()

ang_pix = (1./cdelt2)/60./60.
if margin_unit == 'arcsec':
    margin  = int(margin*ang_pix) # convert from arcsec to pix

# np.abs to make sure it's positive

x_scale   = np.abs(cdelt1*60*60)  # plate scale (dx) [arcsec per pixel]
y_scale   = np.abs(cdelt2*60*60)  # plate scale (dy) [arcsec per pixel]

# read list with targets 

input_file = input_catalogue
hdu        = fits.open(input_file)
header     = hdu[1].header
data       = hdu[1].data
cols= hdu[1].columns
cols_names = np.asarray([cn.lower() for cn in cols.names])

# column names of the input fits table containing the positions
col_name_ra   = cols.names[np.where(cols_names == 'ra')[0][0]]
col_name_dec  = cols.names[np.where(cols_names == 'dec')[0][0]]
if 'unique_id' in cols_names:
    col_name_id   = cols.names[np.where(cols_names == 'unique_id')[0][0]]
else:
    col_name_id   = cols.names[np.where(cols_names == 'id')[0][0]]

col_name_conf = cols.names[np.where(cols_names == 'confidence')[0][0]]

if 'short' in cols_names:
    col_name_class = cols.names[np.where(cols_names == 'short')[0][0]] # column name of the identification of the object
else:
    col_name_class = cols.names[np.where(cols_names == 'lead_line')[0][0]] # column name of the identification of the object

if 'redshift' in cols_names:
    col_name_red = cols.names[np.where(cols_names == 'redshift')[0][0]]
else:
    col_name_red = cols.names[np.where(cols_names == 'z')[0][0]]

if 'lambda_sn' in cols_names:
    col_name_lam = cols.names[np.where(cols_names == 'lambda_sn')[0][0]]
else:
    col_name_lam = cols.names[np.where(cols_names == 'sn')[0][0]]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
try:
    field_numbers = int(field[-2:]) # if field nunber is given in field entry use that
    fieldprovided = True
except:
    field_numbers = (np.asarray([int(str(a)[1:3]) for a in data[col_name_id]])) # otherwise get field number from (MUSE) ids
    fieldprovided = False

if 'candels' in area and 'cdfs' in field:
    here_field = np.where(int(field[5:]) == field_numbers)[0]
    if fieldprovided:
        here_field = np.arange(len(data[col_name_ra]))
elif 'udf' in field:
    here_field = np.where(int(field[-1]) == field_numbers)[0]
elif 'califa' in field:
    here_field = np.arange(len(data[col_name_ra]))
elif 'candels' in area and 'cosmos' in field:
    here_field = np.where(int(field[7:]) == field_numbers)[0]
    if fieldprovided:
        here_field = np.arange(len(data[col_name_ra]))
elif 'cosmos' in area and 'group' in field:
    here_field = np.arange(len(data[col_name_ra]))
else:
    sys.exit(' ERROR: The field name and area combinatio provided (field,area) = ('+field+','+area+') has not "here_field" definitions assigned to them in /Users/kschmidt/work/GitHub/python/run_script_galfit_wrapper.py - add one or change setup file ')

ra_sn = data[col_name_ra][here_field][from_1:to_1]
dec_sn = data[col_name_dec][here_field][from_1:to_1]
obj_ID = np.asarray(data[col_name_id][here_field][from_1:to_1],dtype=str)

try:
    obj_red = data[col_name_red][here_field][from_1:to_1]
except:
    try:
        obj_red = data['Z'][here_field][from_1:to_1]
    except:
        obj_red = data['LAMBDA_SN'][here_field][from_1:to_1]/1215.67-1
try:
    obj_conf = data[col_name_conf][here_field][from_1:to_1]
except:
    obj_conf = np.zeros(len(obj_ID))

try:
    obj_lam = data[col_name_lam][here_field][from_1:to_1]
except:
    obj_lam = (obj_red+1)*1215.67

try:
    obj_class = data[col_name_class][first:last]
except: # if not identified yet and you want to fit all objects
    ident = False

if ident == False:
    obj_class = [False for i in ra_sn]

ra_sn_obj = []
dec_sn_obj = []
ID_obj = []

ra_sn_all_obj = []
dec_sn_all_obj = []
ID_all_obj = [] 

obj_ID_here_old = -1
obj_ID_here_all_old = -1
for i in xrange(len(obj_ID)):
    obj_ID_here_new = obj_ID[i]
    obj_ID_here_all_new = obj_ID[i]
    # list with all objects
    if obj_ID_here_new != obj_ID_here_old:
        ra_sn_all_obj.append(ra_sn[i])
        dec_sn_all_obj.append(dec_sn[i])
        ID_all_obj.append(obj_ID_here_new)
        obj_ID_here_all_old = obj_ID_here_all_new
    # Just fit the objects you want
    if obj_ID_here_new != obj_ID_here_old and obj_class[i]==ident:
        ra_sn_obj.append(ra_sn[i])
        dec_sn_obj.append(dec_sn[i])
        ID_obj.append(obj_ID_here_new)
        obj_ID_here_old = obj_ID_here_new

# star for PSF
if use_for == 'psf':
    ra_sn_obj = [ra_sn_psf]
    dec_sn_obj = [dec_sn_psf]
    ID_obj = [psf_id]   
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# read Guo catalogue to get all visible objects

if Guo_catalogue is not None:
    Guo_open     = open(Guo_catalogue,'r')
    Guo_data     = np.genfromtxt(Guo_open)
    Guo_ID       = [row[0] for row in Guo_data] # ID of Object in cataloge
    Guo_RA       = [row[2] for row in Guo_data]
    Guo_DEC      = [row[3] for row in Guo_data]
    Guo_open.close()
else:
    Guo_ID       = [9999]
    Guo_RA       = [0.0]
    Guo_DEC      = [0.0]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# positions always in pixels for GALFIT!
wcs_obj = wcs.WCS(header_hst, relax=False)

coords = wcs_obj.wcs_world2pix(ra_sn_obj,dec_sn_obj,0)
x_pos = coords[0]
y_pos = coords[1]

if use_for == 'psf':
    try:
        x_pos = [float(x_pos_psf)]
        y_pos = [float(y_pos_psf)]
    except:
        print('   There are no positions for this PSF star yet!')

# positions of all lines (not objects)
coords_all = wcs_obj.wcs_world2pix(ra_sn_all_obj,dec_sn_all_obj,0)
x_pos_all = coords_all[0]
y_pos_all = coords_all[1]

coords_Guo = wcs_obj.wcs_world2pix(Guo_RA,Guo_DEC,0)
x_pos_Guo = coords_Guo[0]
y_pos_Guo = coords_Guo[1]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# find objects close to Lya emitters in Guo catalogue
close_here = [] # array with xy tuples of positions of close objects
close_here_ID = [] # closest object, so hopefully the Guo counterpart
close_here_x = []
close_here_y = []
for i in xrange(len(ID_obj)):
    dist = np.sqrt((x_pos[i]-x_pos_Guo)**2+(y_pos[i]-y_pos_Guo)**2)
    close_bool = dist<margin
    try:
        closest = np.where(dist == np.min(dist))[0][0]
        close_here_ID.append(int(Guo_ID[closest]))
        close_here.append(zip(x_pos_Guo[close_bool],y_pos_Guo[close_bool]))
        close_here_x.append(x_pos_Guo[close_bool])
        close_here_y.append(y_pos_Guo[close_bool])
    except:
        close_here_ID.append([0])
        close_here.append([0,0])
        close_here_x.append([0])
        close_here_y.append([0])

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Read counterparts from Tanya

if counterparts_catalogue is not None:
    input_file_cpart   = counterparts_catalogue
    hdu_cpart          = fits.open(input_file_cpart)
    header_cpart       = hdu_cpart[1].header
    data_cpart         = hdu_cpart[1].data
    cols_cpart         = hdu_cpart[1].columns
    ID_cpart           = data_cpart['UNIQUE_ID']
    Guo_ID_cpart       = data_cpart['GUO_ID']
    Guo_sep_cpart      = data_cpart['GUO_SEP']
    Skelton_ID_cpart   = data_cpart['SKELTON_ID']
    Skelton_sep_scpart = data_cpart['SKELTON_SEP']

    cpart_dict = dict()
    for i in xrange(len(ID_cpart)):
        cpart_dict[ID_cpart[i]] = np.asarray([(Guo_ID_cpart[i]),
                                              float(Guo_sep_cpart[i]),
                                              (Skelton_ID_cpart[i]),
                                              float(Skelton_sep_scpart[i])])
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Read comments
if file_comments!=None:
    
    file_comments_open = open(file_comments,'r')
    file_comments_data = np.genfromtxt(file_comments_open,delimiter='\t',
                                       dtype='str')
    file_comments_ids = np.asarray([int(a[0]) for a in file_comments_data])
    file_comments_comments = np.asarray([a[1] for a in file_comments_data])

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Read parameters from 814 if needed
if use_shape_params != None:
    hdu_814 = fits.open(use_shape_params)
    data_814 = hdu_814[1].data 
    unique_id_814 = np.asarray(data_814['UNIQUE_ID'])
    here_no_multi = np.asarray([i for i in xrange(len(unique_id_814)) 
                     if '_' not in unique_id_814[i]])
    axis_ratio_814 = data_814['axis_ratio'][here_no_multi]
    sersic_exp_814 = data_814['sersic_exp'][here_no_multi]
    magnitude_814 = data_814['magnitude'][here_no_multi]
    comment_814 = data_814['comment'][here_no_multi]
    match_814 = data_814['MATCH'][here_no_multi]
    clumpy_814 = data_814['CLUMPY'][here_no_multi]
    R_e_814 = data_814['R_e'][here_no_multi]
    # convert radius to pixel in this band
    pix_814 = 8.333333e-06
    R_e_814 = R_e_814*(pix_814/cdelt2)
    position_angle_814 = data_814['position_angle'][here_no_multi]
    RA_814 = data_814['RA'][here_no_multi]
    DEC_814 = data_814['DEC'][here_no_multi]
    hdu_814.close()
    coords_814 = wcs_obj.wcs_world2pix(RA_814,DEC_814,0)
    x_814 = coords_814[0]
    y_814 = coords_814[1]

    x_pos = x_814
    y_pos = y_814
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Image region to fit (xmin xmax ymin ymax)
    
xmin = int(np.min(x_pos)-margin)
xmax = int(np.max(x_pos)+margin)
ymin = int(np.min(y_pos)-margin)
ymax = int(np.max(y_pos)+margin)

# get input for function

function = functions()
if obj_type == 'sersic':
    stuff = function.sersic(x_pos)
if obj_type == 'expdisk':
    stuff = function.exponential(x_pos)
if obj_type == 'gaussian':
    stuff = function.gaussian(x_pos)
if obj_type == 'moffat':
    stuff = function.moffat(x_pos)
if obj_type == 'nuker':
    stuff = function.nuker(x_pos)
if obj_type == 'psf':
    stuff = function.psf(x_pos)

fix_other_params = {}
for st in stuff[0].keys():
    fix_other_params[st] = [1 for i in x_pos]

orig_stuff = copy.deepcopy(stuff[0])
params_all_lines = copy.deepcopy(stuff[0])
params_all_lines_errors = copy.deepcopy(stuff[0])
numbers_for_file = stuff[1]

if use_shape_params != None:
    params_all_lines['axis ratio'] = axis_ratio_814
    params_all_lines['sersic exp'] = sersic_exp_814
    params_all_lines['magnitude'] = magnitude_814
    params_all_lines['R_e'] = R_e_814
    params_all_lines['position angle'] = position_angle_814

    fix_other_params['axis ratio'] = [0 for i in x_pos]
    fix_other_params['sersic exp'] = [0 for i in x_pos]
    fix_other_params['magnitude'] = [1 for i in x_pos]
    fix_other_params['R_e'] = [0 for i in x_pos]
    fix_other_params['position angle'] = [0 for i in x_pos]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# copy so params_all_lines is not overwritten
failed_because = copy.deepcopy(stuff[0]) 
for p in failed_because.keys():
    for thing in xrange(len(failed_because[p])):
        failed_because[p][thing] = 0
failed_because_all_lines = failed_because

# margin
margin_all_objs = [margin for i in x_pos]

# array for fixing positions
fix_pos_bool = [fix_pos for i in x_pos]
if fix_pos == False:
    fix_posx = [1 for i in x_pos]
    fix_posy = [1 for i in y_pos]
else:
    fix_posx = [0 for i in x_pos]
    fix_posy = [0 for i in y_pos]

# there seems to be a shift in the positions, so this is necessary
if use_shape_params != None:
    fix_posx = [1 for i in x_pos]
    fix_posy = [1 for i in y_pos]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# parameters for close objects

params_close_here = []
failed_because_close_here = []
fix_other_params_close_here = []
fix_posx_close_here = []
fix_posy_close_here = []

if use_for != 'psf':
    for a in close_here_x:
        if obj_type == 'sersic':
            stuff = function.sersic(a)
        if obj_type == 'expdisk':
            stuff = function.exponential(a)
        params_close_here.append(stuff[0])
        failed_because = copy.deepcopy(stuff[0])
        for p in failed_because.keys():
            for thing in xrange(len(failed_because[p])):
                failed_because[p][thing] = 0
        failed_because_close_here.append(failed_because)
    
        p_fix = {}
        for p in stuff[0].keys():
            p_fix[p] = np.asarray([1 for i in stuff[0]['magnitude']])
        fix_other_params_close_here.append(p_fix)

        fix_posx_close_here.append(np.asarray([fix_posx[0] for i in stuff[0]
                                               ['magnitude']]))
        fix_posy_close_here.append(np.asarray([fix_posy[0] for i in stuff[0]
                                               ['magnitude']]))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# get array of additional objects that were fitted in 814
if use_shape_params != None:
    fit_other_objects = [[] for i in here_no_multi]
    inc_other_objects = [[] for i in here_no_multi]
    per_obj = []
    per_obj_a = []
    i = 0
    for a in xrange(len(unique_id_814)):
        if '_' in unique_id_814[a]:
            per_obj.append(True)
            per_obj_a.append(a)
        if '_' not in unique_id_814[a] or '_' in unique_id_814[a] and \
             i == len(here_no_multi):
            here_multi_ra  = np.asarray(data_814['RA'][per_obj_a])
            here_multi_dec = np.asarray(data_814['DEC'][per_obj_a])
            if len(per_obj_a)>0:
                here_multi_xy = wcs_obj.wcs_world2pix(here_multi_ra,
                                                      here_multi_dec,0)
                close_here_x[i-1] = here_multi_xy[0]
                close_here_y[i-1] = here_multi_xy[1]
                fix_posx_close_here[i-1] = [1 for j in per_obj_a]
                fix_posy_close_here[i-1] = [1 for j in per_obj_a]
            else:
                close_here_x[i-1] = []
                close_here_y[i-1] = []
                fix_posx_close_here[i-1] = []
                fix_posy_close_here[i-1] = []

            fit_other_objects[i-1]  = per_obj
            inc_other_objects[i-1]  = copy.copy(per_obj)
            params_close_here[i-1]['axis ratio'] = data_814['axis_ratio']\
                                                 [per_obj_a]
            params_close_here[i-1]['sersic exp'] = data_814['sersic_exp']\
                                                 [per_obj_a]
            params_close_here[i-1]['R_e'] = data_814['R_e'][per_obj_a]
            pix_814 = 8.333333e-06
            params_close_here[i-1]['R_e'] = params_close_here[i-1]['R_e']*\
                                          (pix_814/cdelt2)
            params_close_here[i-1]['position angle'] = data_814['position_angle'][per_obj_a]
            params_close_here[i-1]['magnitude'] = data_814['magnitude']\
                                                 [per_obj_a]


            failed_because_close_here[i-1]['axis ratio'] = [0 for j in 
                                                            per_obj_a]
            failed_because_close_here[i-1]['sersic exp'] = [0 for j in 
                                                            per_obj_a]
            failed_because_close_here[i-1]['R_e'] = [0 for j in per_obj_a]
            failed_because_close_here[i-1]['magnitude'] = [0 for j in 
                                                           per_obj_a]
            failed_because_close_here[i-1]['position angle'] = [0 for j in 
                                                                per_obj_a]



            fix_other_params_close_here[i-1]['axis ratio'] = [0 for j in 
                                                              per_obj_a]
            fix_other_params_close_here[i-1]['sersic exp'] = [0 for j in 
                                                              per_obj_a]
            fix_other_params_close_here[i-1]['R_e'] = [0 for j in 
                                                       per_obj_a]
            fix_other_params_close_here[i-1]['position angle'] = [0 for j in 
                                                                  per_obj_a]
            fix_other_params_close_here[i-1]['magnitude'] = [1 for j in 
                                                             per_obj_a]

            i+=1
            per_obj = []
            per_obj_a = []

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if use_for != 'psf':
    # array for individual objects in the field of view
    if use_shape_params == None:
        fit_other_objects = [[fit_others for j in close_here_x[i]] 
                             for i in xrange(len(close_here_x))]
        inc_other_objects = [[fit_others for j in close_here_x[i]] 
                             for i in xrange(len(close_here_x))]
else:
    close_here_x = [[x_pos_psf]]
    close_here_y = [[y_pos_psf]]
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def replace(string):
    # (for fixed and failed parameters)
    string = string.replace("[","")
    string = string.replace("]","")
    string = string.replace("*","")
    string = string.replace(",","")
    string = string.replace("(","")
    string = string.replace(")","")
    string = string.replace(" ","")
    return string
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# write information into input parameter files for GALFIT
def write_input_file(params_all_lines):

    for obj in xrange(len(x_pos)):
        xmin = int(float(x_pos[obj])-float(margin_all_objs[obj]))
        xmax = int(float(x_pos[obj])+float(margin_all_objs[obj]))
        ymin = int(float(y_pos[obj])-float(margin_all_objs[obj]))
        ymax = int(float(y_pos[obj])+float(margin_all_objs[obj]))
        
        params_open = open(input_params+'_'+str(ID_obj[obj])+'.input','w')
        
        # This needs to be in the parameter file header
        params_open.write('================================================================================\n')
        params_open.write('# IMAGE and GALFIT CONTROL PARAMETERS\n')
        params_open.write('A) '+input_image+'\n')
        params_open.write('B) '+output_image+'_'+str(ID_obj[obj])+'.fits'+'\n')
        params_open.write('C) '+str(sigma_image)+'\n')
        params_open.write('D) '+str(psf_image)+'\n')
        params_open.write('E) '+str(psf_sampling)+'\n')
        params_open.write('F) '+str(bad_pixel)+'\n')
        params_open.write('G) '+str(param_constr)+'\n')
        params_open.write('H) '+str(xmin)+' '+str(xmax)+' '+str(ymin)+' '+
                          str(ymax)+'\n')
        params_open.write('I) '+str(int(xconv))+' '+str(int(yconv))+'\n')
        params_open.write('J) '+str(phot_zero_ab)+'\n')
        params_open.write('K) '+str(x_scale)+' '+str(y_scale)+'\n')
        params_open.write('O) '+disp_type+'\n')
        params_open.write('P) '+str(options)+'\n\n')
    
        # enter information for objects

        params_open.write('0) '+obj_type+'\n') # Object type
        fix_x = str(fix_posx[obj])
        fix_y = str(fix_posy[obj])

        params_open.write('1) '+str(x_pos[obj])+' '+str(y_pos[obj])+
                          ' '+fix_x+' '+fix_y+' \n')

        params_open.write('Z) 0 \n\n')

        # write parameters for specific function
        for p in xrange(len(params_all_lines.keys())):
            p_here = params_all_lines.keys()[p]
            p_here_value = params_all_lines[p_here][obj]
            p_here_number = numbers_for_file[p_here][0]
            p_here_fix = fix_other_params[p_here][obj]
            params_open.write(str(p_here_number)+') '+str(p_here_value)+
                              ' '+str(p_here_fix)+'\n')
        for other_obj in xrange(len(close_here_x[obj])):
            try:
                if fit_other_objects[obj][other_obj] == True:
                    params_open.write('0) '+obj_type+'\n') # Object type
                    
                    fix_x_close_here = str(int(fix_posx_close_here[obj][other_obj]))
                    fix_y_close_here = str(int(fix_posy_close_here[obj][other_obj]))

                    write_fix = ' '+fix_x_close_here+' '+fix_y_close_here+' \n'

                    params_open.write('1) '+str(close_here_x[obj][other_obj])+
                                      ' '+str(close_here_y[obj][other_obj])+
                                      write_fix)

                    for p in xrange(len(params_all_lines.keys())):
                        p_here = params_close_here[obj].keys()[p]
                        p_here_value = params_close_here[obj][p_here][other_obj]
                        p_here_number = numbers_for_file[p_here][0]
                        p_fix = fix_other_params_close_here[obj][p_here][other_obj]
                        params_open.write(str(p_here_number)+') '+
                                          str(p_here_value)+' '+str(p_fix)+'\n')
                    params_open.write('Z) 0 \n\n')
            except:
                pass

        if fix_sky == True:
            fixing_sky = 0
        else:
            fixing_sky = 1
        if sky_value != None:
            params_open.write('# sky\n\n') # Sky
            params_open.write('0) sky\n')
            params_open.write('1) '+str(sky_value)+' '+str(fixing_sky)+' \n')
            #params_open.write('2) '+str(sky_dx)+' '+str(fixing_sky)+' \n')
            #params_open.write('3) '+str(sky_dy)+' '+str(fixing_sky)+' \n')
            # no variation in direction
            params_open.write('2) '+str(sky_dx)+' '+str(0)+' \n')
            params_open.write('3) '+str(sky_dy)+' '+str(0)+' \n')
            params_open.write('Z) 0 \n\n')

        params_open.close()

write_input_file(params_all_lines)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### PART 2 - run GALFIT

failed = []
failed_bool = [False for i in x_pos]
bombed_bool = [False for i in x_pos]

try:
    match_clear = copy.copy(match_814) 
    match_clumpy = [False for i in x_pos]
    here_clumpy_814 = np.where(clumpy_814==1)[0]
    for i in here_clumpy_814:
        match_clumpy[i] = True
except:
    match_clear = ['multi' for i in x_pos] # set default tick clumpy/clear/nothing/multiple
    match_clumpy = [False for i in x_pos]

own_comments = ['-' for i in x_pos] 
new_objects_x  = []
new_objects_y  = []
new_objects_ids = []
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def run_galfit(redo_bool,obj):

    print('   Now looking at '+str(len(x_pos))+' objects...')

    all_lines = []
    all_lines_others = [] # for other objects

    # reset failed and bombed for each run
    failed_bool[obj] = False
    bombed_bool[obj] = False
        
    if redo_bool[obj]==True:
        # TODO: Why do I need the full path here? Can I change that?
        command = (galfit_path+' -outsig '+input_params+'_'+
                   str(ID_obj[obj])+'.input')  
        print('   '+str(command))
        os.system(command)

        # Read and sort the output files
        
        try:
            log_open = open('fit.log','r')
            lines = log_open.readlines()
            # since there are two * in the log file if the fit went wrong
            first = True 
            other_lines = []

            for line in lines:
                if obj_type in line and first == False:
                    other_lines.append(line)
                if obj_type in line and first == True:
                    all_lines.append(line)
                    first = False # the first line is for the main object
            all_lines_others.append(other_lines)
            log_open.close()
            if save_stuff == True:
                shutil.copy2('fit.log',out_log_dir+'/fit_'+band+'_'+str(ID_obj[obj])+
                             '.log')
            # fit.log is always appended and therefore hard to machine-read
            os.remove('fit.log')
        except IOError:
            bombed_bool[obj] = True
            failed_bool[list_j[obj]] = True
            all_lines.append('bombed')
            all_lines_others.append(['bombed'])
    else: # if not to be redone, append empty dummy list
        all_lines_others.append([])
        all_lines.append([])

    # read the last parameters of the lines where it worked

    for p in params_all_lines.keys():
        position_here = numbers_for_file[p][1]
        p_array = []
        
        # for main object
        if bombed_bool[obj]==False and redo_bool[obj]==True:
            line = all_lines[0].split()
            value_here = line[position_here]

            params_all_lines[p][list_j[obj]] = float(replace(value_here))

            if '*' in value_here:
                failed_bool[list_j[obj]] = True
                failed_because_all_lines[p][list_j[obj]] = 1.
            else:
                failed_because_all_lines[p][list_j[obj]] = 0.

            if use_for != 'psf':
                # for other objects
                l = 0
                for fit_o in xrange(len(fit_other_objects[obj])):
                    if fit_other_objects[obj][fit_o] == True:
                        line = all_lines_others[0][l].split()
                        value_here = line[position_here]
                        val = float(replace(value_here))
                        if  '*' in value_here:
                            failed_because_close_here[list_j[obj]][p][fit_o]=1.
                        else:
                            failed_because_close_here[list_j[obj]][p][fit_o]=0.
                        params_close_here[list_j[obj]][p][fit_o] = val
                        l+=1

    # NOTE: We don't use the new positions, because they are often good enough
    # as a first guess to get them right in one run, but if they were wrong
    # the first run, it's no use to start over with the wrong positions in the
    # second run.

    # move galfit files 
    galfit_files = [f for f in listdir('.') if 'galfit.' in f]
    for f in galfit_files:
        shutil.move(f,out_log_dir+f)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### PART 3 - use a GUI for easier iteration

# values for the cubehelix colormap
ch_g = 1.0
ch_s = 1.7
ch_r = -1.5
ch_h = 0.7

if counterparts_catalogue == None:
    add_width = 0
else:
    add_width = 400

# list for objects that failed but looked ok
failed = []

def cubehelix(gamma=1.0, s=0.5, r=-1.5, h=1.0):
    """
    creates a lookup table for the cubehelix colormap
    """
    def get_color_function(p0, p1):
        def color(x):
            xg = x ** gamma
            a = h * xg * (1 - xg) / 2
            phi = 2 * np.pi * (s / 3 + r * x)
            return xg + a * (p0 * np.cos(phi) + p1 * np.sin(phi))
        return color

    array = np.empty((256, 3))
    abytes = np.arange(0, 1, 1/256.)
    array[:, 0] = get_color_function(-0.14861, 1.78277)(abytes) * 255
    array[:, 1] = get_color_function(-0.29227, -0.90649)(abytes) * 255
    array[:, 2] = get_color_function(1.97294, 0.0)(abytes) * 255
    return array

def vminvmax(image, scale=0.95, nans=0):
    """
    vmin,vmax = vminvmax(data, scale=0.99, nans=0)

    In:
    ---
    data - 2D array (intensity image)
    scale (=0.99) - scale factor for colorbar
    nans (=0) value with wich NaNs in the array are replaced

    Out:
    ---
    vmin, vmax -- values to be used in matplotlib.imshow as vmin & vmax
    """
    assert scale <= 1 and scale > 0 
    data = image.copy()    

    data[np.isnan(data)] = nans
    if len(data.shape) > 1:
        data = data.flatten()
    data = np.sort(data)

    length = len(data)

    vmax_index = int(np.floor(scale*length))
    vmin_index = int(np.ceil((1-scale)*length))

    vmin = data[vmin_index]
    vmax = data[vmax_index]

    return vmin, vmax

redo_bool = [True for i in x_pos] # contains True if redo and False if not
list_j = [i for i in xrange(len(redo_bool)) if redo_bool[i]==True]

run_galfit(redo_bool,0)
j=0

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class main_Window(QtGui.QMainWindow):

    def __init__(self):
        super(main_Window, self).__init__()

        self.initUI()

    def initUI(self):

        self.list_j = [i for i in xrange(len(redo_bool)) if 
                       redo_bool[i]==True]
        self.j = j

        self.statusBar().showMessage('')
        self.main_GUI = main_GUI()
        self.setCentralWidget(self.main_GUI)

        self.main_GUI.signal_close[bool].connect(self.close_now)
    
        self.update(self.j)
        self.setGeometry(800, 300, 850+add_width, 450)
        self.main_GUI.signal_ID.connect(self.update)

        self.center()
        self.show()

    def update(self, new_j):
        self.j = new_j
        ID = ID_obj[self.list_j[self.j]]
        self.setWindowTitle('Looking at object ID '+str(ID)) 

    def center(self):
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()     
        qr.moveCenter(cp)
        self.move(qr.topLeft())

    def close_now(self):
        list_j = [i for i in xrange(len(redo_bool)) if 
                  redo_bool[i]==True]
        write_input_file(params_all_lines)
        if len(list_j) != 0:
            self.j = 0
            # move output images to different folder
            for files in os.listdir(image_dir):
                if 'imgblock' in files:
                    if save_stuff == True:
                        shutil.copy2(image_dir+files,result_folder)
                    os.remove(image_dir+files)
            run_galfit(redo_bool)
            # do again if there is something to redo
            self.initUI()
        else:
            self.close()
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class main_GUI(QtGui.QWidget):
   
    signal_close = QtCore.pyqtSignal(bool)
    signal_ID = QtCore.pyqtSignal(int)

    def __init__(self):
        super(main_GUI, self).__init__()

        self.initUI(j)
        
    def initUI(self,j):

        self.list_j = [i for i in xrange(len(redo_bool)) if 
                       redo_bool[i]==True]
        self.j = j
        self.c = 999

        # set line edits with parameters for the fit
        self.param_dict = dict()
        self.cb_params = dict()
        self.fail_labels = dict()
        for param in xrange(len(params_all_lines.keys())):
            param_here = params_all_lines.keys()[param] 
            self.le_params = QtGui.QLineEdit(self)
            self.le_params.setGeometry(add_width+50+80*param,320,50,25)
            ind_here = self.list_j[self.j]
            self.le_params.setText(str(params_all_lines[param_here][ind_here]))
            params_label = QtGui.QLabel(param_here, self)
            params_label.move(add_width+50+80*param,300)
            self.param_dict[param_here] = self.le_params

            # checkbox for each object for each parameter
            self.cb_fix_each = QtGui.QCheckBox('',self)
            self.cb_fix_each.setToolTip('Fix this parameter. \n axis ratio: min=0.1, fix=1 \n sersic exp: min=0.2, max=10, fix=1 \n magnitude: min=29, fix=29 \n Radius: min=0.5, fix=1')
            self.cb_fix_each.move(add_width+50+80*param,350)
            self.cb_params[param_here] = self.cb_fix_each

            # Information if it failed or bombed for each object
            self.fail_label_all = QtGui.QLabel('                 ', self)
            self.fail_label_all.setGeometry(add_width+50+80*param,370,100,10)
            self.fail_labels[param_here] = self.fail_label_all
        self.write_fail_all()
            
        self.fix_label = QtGui.QLabel('fix', self)
        self.fix_label.move(add_width+15,355)

        # comment line
        own_comm_title = QtGui.QLabel('You can enter a comment:', self)
        own_comm_title.move(add_width,395)
        self.own_comm = QtGui.QLineEdit(self)
        self.own_comm.setGeometry(add_width,410,350,20)
        self.own_comm.setText(own_comments[self.j])

        # checkbox for if the object is clumpy
        self.cb_clumpy = QtGui.QCheckBox('Clumpy',self)
        self.cb_clumpy.setToolTip('This object is clumpy.')
        self.cb_clumpy.move(add_width+650,200)
        self.cb_clumpy.stateChanged.connect(self.write_clumpy)

        # checkbox for if the object is obvious
        self.cb_clear = QtGui.QCheckBox('Clear',self)
        self.cb_clear.setToolTip('This object can be clearly matched.')
        self.cb_clear.move(add_width+650,220)
        self.cb_clear.stateChanged.connect(self.write_clear)

        # checkbox for if the object is not there
        self.cb_nothing = QtGui.QCheckBox('Nothing',self)
        self.cb_nothing.setToolTip('There is nothing there.')
        self.cb_nothing.move(add_width+650,240)
        self.cb_nothing.stateChanged.connect(self.write_nothing)

        # checkbox for if there are multiple objects 
        self.cb_multi = QtGui.QCheckBox('Multiple',self)
        self.cb_multi.setToolTip('There are multiple objects. If this box is checked, the parameters of all objects that are fitted will be saved.')
        self.cb_multi.move(add_width+650,260)
        self.cb_multi.stateChanged.connect(self.write_multi)

        self.toggle_match()
                                                        
        # line edit for the margin
        self.le_margin = QtGui.QLineEdit(self)
        self.le_margin.setGeometry(add_width+650,170,50,25)
        self.le_margin.setText(str(margin_all_objs[self.list_j[self.j]]))
        margin_label = QtGui.QLabel("Margin", self)
        margin_label.move(add_width+710,175)

        long_line = QtGui.QLabel('___________________________', self)
        long_line.move(add_width+650,90)

        # line edit for positions of additional objects
        self.le_addx = QtGui.QLineEdit(self)
        self.le_addx.setGeometry(add_width+650,125,45,25)
        self.le_addx.setText('0')
        add_labelx = QtGui.QLabel("Additional object", self)
        add_labelx.move(add_width+650,105)
        self.le_addy = QtGui.QLineEdit(self)
        self.le_addy.setGeometry(add_width+710,125,45,25)
        self.le_addy.setText('0')

        def add_object():
            h = self.list_j[self.j]
            new_posx = self.le_addx.text()
            new_posy = self.le_addy.text()
            new_objects_x.append(int(new_posx))
            new_objects_y.append(int(new_posy))
            new_objects_ids.append(int(ID_obj[h].split('.')[0]))
            here_oldx = close_here_x[h]
            here_oldx = np.append(here_oldx,int(new_posx))
            close_here_x[h] = here_oldx
            here_oldy = close_here_y[h]
            here_oldy = np.append(here_oldy,int(new_posy))
            close_here_y[h] = here_oldy

            new_here = np.append(fit_other_objects[h],True)
            fit_other_objects[h] = new_here
            inc_other_objects[h] = copy.copy(new_here)
            new_herex = np.append(fix_posx_close_here[h],True)
            fix_posx_close_here[h] = new_herex
            new_herey = np.append(fix_posy_close_here[h],True)
            fix_posy_close_here[h] = new_herey

            for p in xrange(len(params_all_lines.keys())):
                p_here2 = params_close_here[h].keys()[p]
                old_here2 = params_close_here[h][p_here2]
                new_param2 = np.append(old_here2,orig_stuff[p_here2][0])
                params_close_here[h][p_here2] = new_param2
                old_here2 = fix_other_params_close_here[h][p_here2]
                new_here2 = np.append(old_here2,1)
                fix_other_params_close_here[h][p_here2] = new_here2
                
                old_here2 = failed_because_close_here[h][p_here2]
                new_here2 = np.append(old_here2,False)
                failed_because_close_here[h][p_here2] = new_here2

            self.le_addx.setText('0')
            self.le_addy.setText('0')

        button_k = QtGui.QPushButton('+', self)
        button_k.setToolTip("Click to add object with this position to list.")
        #button_k.move(add_width+760,125)
        button_k.setGeometry(add_width+760,120,50,30)
        button_k.clicked.connect(add_object)

        long_line = QtGui.QLabel('___________________________', self)
        long_line.move(add_width+650,145)

        # Information if it failed or bombed
        self.fail_label = QtGui.QLabel('                 ', self)
        self.fail_label.setGeometry(add_width+740,80,100,10)
        self.write_fail()

        # Information for redshift
        self.red_label = QtGui.QLabel('redsift: '+
                                      str(obj_red[self.list_j[self.j]]),self)
        self.red_label.setGeometry(add_width+650,375,100,10)

        # Information for confidence
        self.conf_label = QtGui.QLabel('confidence: '+
                                       str(obj_conf[self.list_j[self.j]]),self)
        self.conf_label.setGeometry(add_width+650,395,100,10)

        # Information for lambda
        self.lam_label = QtGui.QLabel('lambda: '+
                                      str(obj_lam[self.list_j[self.j]]),self)
        self.lam_label.setGeometry(add_width+650,415,100,10)
        self.write_red_conf_lam_label()

        # set line edits with positions
        self.posx = QtGui.QLineEdit(self)
        self.posx.setGeometry(add_width+500,320,70,25)
        self.posx.setText(str(np.round(x_pos[self.list_j[self.j]],2)))
        posx_label = QtGui.QLabel('x pos', self)
        posx_label.move(add_width+500,300)

        self.posy = QtGui.QLineEdit(self)
        self.posy.setGeometry(add_width+570,320,70,25)
        self.posy.setText(str(np.round(y_pos[self.list_j[self.j]],2)))
        posy_label = QtGui.QLabel('y pos', self)
        posy_label.move(add_width+570,300)

        # checkbox for each object for positions
        self.cb_fix_posx_each = QtGui.QCheckBox('',self)
        self.cb_fix_posx_each.setToolTip('Fix x-position.')
        self.cb_fix_posx_each.move(add_width+500,350)
        self.cb_fix_posy_each = QtGui.QCheckBox('',self)
        self.cb_fix_posy_each.setToolTip('Fix y-position.')
        self.cb_fix_posy_each.move(add_width+570,350)

        # check or uncheck checkboxes for first object
        if self.cb_fix_posx_each.isChecked() and \
           fix_posx[self.list_j[self.j]] != 0:
            self.cb_fix_posx_each.toggle()
        elif not self.cb_fix_posx_each.isChecked() and \
             fix_posx[self.list_j[self.j]] == 0:
            self.cb_fix_posx_each.toggle()

        if self.cb_fix_posy_each.isChecked() and \
           fix_posy[self.list_j[self.j]] != 0:
            self.cb_fix_posy_each.toggle()
        elif not self.cb_fix_posy_each.isChecked() and \
             fix_posy[self.list_j[self.j]] == 0:
            self.cb_fix_posy_each.toggle()

        long_line = QtGui.QLabel('__________________________________________________________________________________________________________________________________________________________________________________________________________', self)
        long_line.move(0,270)

        button_redo = QtGui.QPushButton('Redo', self)
        button_redo.setToolTip('Redo the fit for this object, look at next object.')
        button_redo.move(add_width+650,30)
        button_redo.clicked.connect(self.write_params)
        button_redo.clicked.connect(self.redo)
        button_redo.setShortcut('r')

        button_ok = QtGui.QPushButton('Ok', self)
        button_ok.setToolTip("Don't redo the fit for this object, look at next object.")
        button_ok.move(add_width+650,66)
        button_ok.clicked.connect(self.write_params)
        button_ok.clicked.connect(self.ok)
        button_ok.setShortcut('o')

        img1_label = QtGui.QLabel('Original', self)
        img1_label.move(add_width+50,250)
        img2_label = QtGui.QLabel('Model', self)
        img2_label.move(add_width+250,250)
        img3_label = QtGui.QLabel('Subtracted', self)
        img3_label.move(add_width+450,250)

        self.img1 = pg.GraphicsLayoutWidget(self)
        self.img1.setGeometry(add_width+50,50, 180, 180)
        self.img2 = pg.GraphicsLayoutWidget(self)
        self.img2.setGeometry(add_width+250,50, 180, 180)
        self.img3 = pg.GraphicsLayoutWidget(self)
        self.img3.setGeometry(add_width+450,50, 180, 180)

        # image of counterparts
        if counterparts_catalogue != None and use_for != 'psf':
            self.img0 = pg.GraphicsLayoutWidget(self)
            self.img0.setGeometry(40,0, 260, 260)
            img0_label = QtGui.QLabel('Counterpart', self)
            img0_label.move(40,265)

            # Information for IDs for counterparts
            ID = ID_obj[self.list_j[self.j]]
            try:
                self.cpart_label0 = QtGui.QLabel('Guo ID: '+
                                                 str(cpart_dict[int(ID)][0]),
                                                 self)
                self.cpart_label1 = QtGui.QLabel('Guo sep.: '+
                                                 str(cpart_dict[int(ID)][1]),
                                                 self)
                self.cpart_label2 = QtGui.QLabel('Skelton ID: '+
                                                 str(cpart_dict[int(ID)][2]),
                                                 self)
                self.cpart_label3 = QtGui.QLabel('Skelton sep.: '+
                                                 str(cpart_dict[int(ID)][3]),
                                                 self)
            except:
                self.cpart_label0 = QtGui.QLabel('Guo ID: ',self)
                self.cpart_label1 = QtGui.QLabel('Guo sep.: ',self)
                self.cpart_label2 = QtGui.QLabel('Skelton ID: ',self)
                self.cpart_label3 = QtGui.QLabel('Skelton sep.: ',self)

            self.cpart_label0.setGeometry(40,295,100,12)
            self.cpart_label1.setGeometry(40,312,100,12)
            self.cpart_label2.setGeometry(40,329,100,12)
            self.cpart_label3.setGeometry(40,346,100,12)
            


            self.write_cpart_label()

        # comments
        if file_comments != None:
            ID = int(ID_obj[self.list_j[self.j]])
            try:
                here_ID_comm = [i for i in xrange(len(file_comments_ids))
                                if file_comments_ids[i]==ID][0]
                self.comm_label = QtGui.QLabel('Comment: '+
                                               file_comments_comments[
                                                   here_ID_comm],self)
            except:
                self.comm_label = QtGui.QLabel('Comment: ',self)
            self.comm_label.setGeometry(40,367,400,12)
            self.write_comm_label()


        # checkbox for each object
        self.cb_fit_each = QtGui.QCheckBox('Fit object',self)
        self.cb_fit_each.setToolTip('Fit this object.')
        self.cb_fit_each.move(add_width+650,320)
        self.cb_fit_each.stateChanged.connect(self.fit_this)
        
        # checkbox for each object to include
        self.cb_inc_each = QtGui.QCheckBox('Include object',self)
        self.cb_inc_each.setToolTip('Include this object.')
        self.cb_inc_each.move(add_width+650,340)
        self.cb_inc_each.stateChanged.connect(self.inc_this)
        

        self.display_images()

    def write_clumpy(self,state):
        if state == QtCore.Qt.Checked:
            match_clumpy[self.list_j[self.j]] = True
        else:
            match_clumpy[self.list_j[self.j]] = False
    def write_clear(self,state):
        if state == QtCore.Qt.Checked:
            match_clear[self.list_j[self.j]] = 'clear' 
    def write_nothing(self,state):
        if state == QtCore.Qt.Checked:
            match_clear[self.list_j[self.j]] = 'nothing' 
    def write_multi(self,state):
        if state == QtCore.Qt.Checked:
            match_clear[self.list_j[self.j]] = 'multi' 

    def toggle_match(self):

        if not self.cb_clumpy.isChecked() and match_clumpy[self.list_j
                                                           [self.j]] ==True:
            self.cb_clumpy.toggle()
        elif self.cb_clumpy.isChecked():
            self.cb_clumpy.toggle()

        if not self.cb_clear.isChecked() and match_clear[self.list_j
                                                         [self.j]] =='clear':
            self.cb_clear.toggle()
        elif self.cb_clear.isChecked():
            self.cb_clear.toggle()
                
        if not self.cb_nothing.isChecked() and match_clear[self.list_j
                                                           [self.j]]=='nothing':
            self.cb_nothing.toggle()
        elif self.cb_nothing.isChecked():
            self.cb_nothing.toggle()
        if not self.cb_multi.isChecked() and match_clear[self.list_j
                                                         [self.j]] =='multi':
            self.cb_multi.toggle()
        elif self.cb_multi.isChecked():
            self.cb_multi.toggle()


    def write_fail(self):
        self.fail_label.clear()
        if failed_bool[self.list_j[self.j]] == True:
            self.fail_label.setText('failed!')

    def write_red_conf_lam_label(self):
        self.red_label.clear()
        self.red_label.setText('redshift: '+str(np.round(obj_red[self.list_j[self.j]],2)))
        self.conf_label.clear()
        self.conf_label.setText('confidence: '+str(obj_conf[self.list_j
                                                            [self.j]]))
        self.lam_label.clear()
        self.lam_label.setText('lambda: '+str(np.round(obj_lam[self.list_j[self.j]])))


    def write_cpart_label(self):
        if counterparts_catalogue != None and use_for != 'psf':
            self.cpart_label0.clear()
            self.cpart_label1.clear()
            self.cpart_label2.clear()
            self.cpart_label3.clear()

        ID = ID_obj[self.list_j[self.j]]

        if counterparts_catalogue != None and use_for != 'psf':
            try:
                self.cpart_label0.setText('Guo ID: '+str(cpart_dict[str(ID)][0]))
                self.cpart_label1.setText('Guo sep.: '+str(cpart_dict[str(ID)][1]))
                self.cpart_label2.setText('Skelton ID: '+str(cpart_dict[str(ID)][2]))
                self.cpart_label3.setText('Skelton sep.: '+str(cpart_dict[str(ID)][3]))
            except:
                pass



    def write_comm_label(self):
        if file_comments != None:
            self.comm_label.clear()
            ID = int(ID_obj[self.list_j[self.j]])
            try:
                here_ID_comm = [i for i in xrange(len(file_comments_ids)) if 
                                file_comments_ids[i]==ID][0]
                self.comm_label.setText('Comment: '+
                                        file_comments_comments[here_ID_comm])
            except:
                self.comm_label.setText('Comment: ')


    def write_fail_all(self):
        
        h = self.list_j[self.j]

        for param in xrange(len(params_all_lines.keys())):
            param_here = params_all_lines.keys()[param] 
            self.fail_labels[param_here].clear()

            if self.c == 999:
                if failed_because_all_lines[param_here][h]==1:
                    self.fail_labels[param_here].setText('failed!')
            else:
                if failed_because_close_here[h][param_here][self.c]==1:
                    self.fail_labels[param_here].setText('failed!')

    def fit_this(self,state):
        h = self.list_j[self.j]
        if self.c != 999:
            if state == QtCore.Qt.Checked:
                fit_other_objects[h][self.c] = True
            else:
                fit_other_objects[h][self.c] = False
                # don't include object if it's not fitted
                self.cb_inc_each.setChecked(False)

    def inc_this(self,state):
        h = self.list_j[self.j]
        if self.c != 999:
            if state == QtCore.Qt.Checked:
                inc_other_objects[h][self.c] = True
            else:
                inc_other_objects[h][self.c] = False

    def write_params(self):
        
        h = self.list_j[self.j]

        if self.c == 999:
            for param in xrange(len(params_all_lines.keys())):
                param_here = params_all_lines.keys()[param] 
                line_text = self.param_dict[param_here].text()
                params_all_lines[param_here][h] = str(line_text)

                if self.cb_params[param_here].isChecked():
                    fix_other_params[param_here][h] = 0
                else:
                    fix_other_params[param_here][h] = 1

                # fix or unfix position
                if self.cb_fix_posx_each.isChecked():
                    fix_posx[h] = 0
                else:
                    fix_posx[h] = 1
                if self.cb_fix_posy_each.isChecked():
                    fix_posy[h] = 0
                else:
                    fix_posy[h] = 1

            # write positions
            x_pos[h] = self.posx.text()
            y_pos[h] = self.posy.text()

        else:
            for param in xrange(len(params_all_lines.keys())):
                param_here = params_all_lines.keys()[param] 
                line_text = self.param_dict[param_here].text()
                params_close_here[h][param_here][self.c] = str(line_text)

                if self.cb_params[param_here].isChecked():
                    fix_other_params_close_here[h][param_here][self.c] = 0
                else:
                    fix_other_params_close_here[h][param_here][self.c] = 1

                # fix or unfix positions for other objects
                if self.cb_fix_posx_each.isChecked():
                    fix_posx_close_here[h][self.c] = 0
                else:
                    fix_posx_close_here[h][self.c] = 1
                if self.cb_fix_posy_each.isChecked():
                    fix_posy_close_here[h][self.c] = 0
                else:
                    fix_posy_close_here[h][self.c] = 1

            # write positions
            close_here_x[h][self.c] = self.posx.text()
            close_here_y[h][self.c] = self.posy.text()

    def ok(self):

        self.c = 999
        redo_bool[self.list_j[self.j]] = False
        if failed_bool[self.list_j[self.j]] == True:
            failed.append(ID_obj[self.list_j[self.j]])

        # save comment
        comment = self.own_comm.text()
        own_comments[self.j] = str(comment)

        self.j +=1
        
        if self.j >= len(self.list_j):
            self.signal_close.emit(True)
        else:  
            self.own_comm.setText(own_comments[self.j])
            # check or uncheck checkboxes for next object
            if self.cb_fix_posx_each.isChecked() and \
               fix_posx[self.list_j[self.j]] != 0:
                self.cb_fix_posx_each.toggle()
            elif not self.cb_fix_posx_each.isChecked() and \
                 fix_posx[self.list_j[self.j]] == 0:
                self.cb_fix_posx_each.toggle()

            if self.cb_fix_posy_each.isChecked() and \
               fix_posy[self.list_j[self.j]] != 0:
                self.cb_fix_posy_each.toggle()
            elif not self.cb_fix_posy_each.isChecked() and \
                 fix_posy[self.list_j[self.j]] == 0:
                self.cb_fix_posy_each.toggle()

            write_input_file(params_all_lines)
            run_galfit(redo_bool,self.list_j[self.j])

            if not self.cb_clumpy.isChecked() and match_clumpy[self.list_j
                                                               [self.j]] ==True:
                self.cb_clumpy.toggle()
            elif self.cb_clumpy.isChecked() and match_clumpy[self.list_j
                                                             [self.j]] !=True:
                self.cb_clumpy.toggle()

            if not self.cb_clear.isChecked() and match_clear[self.list_j
                                                             [self.j]]=='clear':
                self.cb_clear.toggle()
            elif self.cb_clear.isChecked() and match_clear[self.list_j
                                                           [self.j]] !='clear':
                self.cb_clear.toggle()
                
            if not self.cb_nothing.isChecked() and \
               match_clear[self.list_j[self.j]]=='nothing':
                self.cb_nothing.toggle()
            elif self.cb_nothing.isChecked() and \
                 match_clear[self.list_j[self.j]]!='nothing':
                self.cb_nothing.toggle()
            if not self.cb_multi.isChecked() and \
               match_clear[self.list_j[self.j]] =='multi':
                self.cb_multi.toggle()
            elif self.cb_multi.isChecked() and \
                 match_clear[self.list_j[self.j]] !='multi':
                self.cb_multi.toggle()

            self.posx.setText(str(np.round(x_pos[self.list_j[self.j]],2)))
            self.posy.setText(str(np.round(y_pos[self.list_j[self.j]],2)))
            self.le_margin.setText(str(margin_all_objs[self.list_j[self.j]]))
            self.signal_ID.emit(self.j)
            self.display_images()
            self.write_fail()
            self.write_fail_all()
            self.write_red_conf_lam_label()
            if use_for != 'psf':
                self.write_cpart_label()
                self.write_comm_label()


    def redo(self):

        write_input_file(params_all_lines)
        run_galfit(redo_bool,self.list_j[self.j])

        self.write_fail()
        self.write_fail_all()
        self.write_red_conf_lam_label()
        if use_for != 'psf':
            self.write_cpart_label()
            self.write_comm_label()
        self.c = 999
        redo_bool[self.list_j[self.j]] = True
        margin_all_objs[self.list_j[self.j]] = self.le_margin.text()

        #self.j +=1

        # uncheck checkboxes for next object

        if self.cb_fix_posx_each.isChecked() and \
           fix_posx[self.list_j[self.j]] != 0:
            self.cb_fix_posx_each.toggle()
        elif not self.cb_fix_posx_each.isChecked() and \
             fix_posx[self.list_j[self.j]] == 0:
            self.cb_fix_posx_each.toggle()

        if self.cb_fix_posy_each.isChecked() and \
           fix_posy[self.list_j[self.j]] != 0:
            self.cb_fix_posy_each.toggle()
        elif not self.cb_fix_posy_each.isChecked() and \
             fix_posy[self.list_j[self.j]] == 0:
            self.cb_fix_posy_each.toggle()

        if self.j >= len(self.list_j):
            self.signal_close.emit(True)
        else:  

            if not self.cb_clumpy.isChecked() and match_clumpy[self.list_j
                                                               [self.j]] ==True:
                self.cb_clumpy.toggle()
            elif self.cb_clumpy.isChecked() and match_clumpy[self.list_j
                                                             [self.j]] !=True:
                self.cb_clumpy.toggle()

            if not self.cb_clear.isChecked() and \
               match_clear[self.list_j[self.j]] =='clear':
                self.cb_clear.toggle()
            elif self.cb_clear.isChecked() and \
                 match_clear[self.list_j[self.j]] !='clear':
                self.cb_clear.toggle()
                
            if not self.cb_nothing.isChecked() and \
               match_clear[self.list_j[self.j]]=='nothing':
                self.cb_nothing.toggle()
            elif self.cb_nothing.isChecked() and \
                 match_clear[self.list_j[self.j]]!='nothing':
                self.cb_nothing.toggle()
            if not self.cb_multi.isChecked() and \
               match_clear[self.list_j[self.j]] =='multi':
                self.cb_multi.toggle()
            elif self.cb_multi.isChecked() and \
                 match_clear[self.list_j[self.j]] !='multi':
                self.cb_multi.toggle()

            self.posx.setText(str(np.round(float(x_pos[self.list_j[self.j]]),2)))
            self.posy.setText(str(np.round(float(y_pos[self.list_j[self.j]]),2)))
            self.le_margin.setText(str(margin_all_objs[self.list_j[self.j]]))
            self.signal_ID.emit(self.j)
            self.display_images()
            self.write_fail()
            self.write_fail_all()
            self.write_red_conf_lam_label()
            if use_for != 'psf':
                self.write_cpart_label()
                self.write_comm_label()

    def display_images(self):

        id_here = ID_obj[self.list_j[self.j]]

        for param in xrange(len(params_all_lines.keys())):
            param_here = params_all_lines.keys()[param] 
            self.param_dict[param_here].setText(str(params_all_lines[param_here]
                                                    [self.list_j[self.j]]))

            # reset fix ceck boxes
            this_param_fix=fix_other_params[param_here][self.list_j[self.j]]
            if self.cb_params[param_here].isChecked() and \
               this_param_fix==1:
                self.cb_params[param_here].toggle()
            elif not self.cb_params[param_here].isChecked() and \
                 this_param_fix==0:
                self.cb_params[param_here].toggle()


        h = self.list_j[self.j]
    
        self.cb_fit_each.blockSignals(True)
        self.cb_inc_each.blockSignals(True)
        self.cb_fit_each.setChecked(True)
        self.cb_inc_each.setChecked(True)
        self.cb_fit_each.blockSignals(False)
        self.cb_inc_each.blockSignals(False)
        
        ## Create image item

        try:
            image_open  = fits.open(output_image+'_'+str(id_here)+".fits")
            image_data1 = image_open[1].data
            image_data2 = image_open[2].data
            image_data3 = image_open[3].data
        except: # in case it bombed there's no fits file
            image_open  = fits.open(placeholder_bombed)
            image_data1 = image_open[0].data
            image_data2 = image_open[0].data
            image_data3 = image_open[0].data

        image_open.close()

        self.img1.clear()
        self.img2.clear()
        self.img3.clear()

        if counterparts_catalogue != None and use_for != 'psf':
            self.img0.clear()

        # image for counterparts
        if counterparts_catalogue != None and use_for != 'psf':
            ID = ID_obj[self.list_j[self.j]]
            if 'III' in counterparts_suff:
                image_name0 = counterparts_imgs+counterparts_suff.replace('III',str(ID))
            else:
                image_name0 = counterparts_imgs+str(ID)+counterparts_suff
            image_data0 = np.array(mpimg.imread(image_name0),dtype='float64')
            img_0 = pg.ImageItem()
            img_0.setImage(np.rot90(np.rot90(np.rot90(image_data0))))
            self.view_img0 = self.img0.addViewBox()
            self.view_img0.addItem(img_0)
            if 'III' in counterparts_suff:
                pass
                #self.view_img0.setRange(QtCore.QRectF(0,0,500,500))
            else:
                self.view_img0.setRange(QtCore.QRectF(300,300,400,400))

        level_min, level_max = vminvmax(image_data1)

        margin_here = margin_all_objs[self.list_j[self.j]]

        img_1 = pg.ImageItem()
        img_1.setImage(np.rot90(image_data1)[::-1])
        img_1.setLookupTable(cubehelix(ch_g,ch_s,ch_r,ch_h))
        img_1.setLevels([level_min,level_max])
        self.view_img1 = self.img1.addViewBox()
        self.view_img1.addItem(img_1)
        pixwidth = int(float(margin_here)*2)
        self.view_img1.setRange(QtCore.QRectF(0,0,pixwidth,pixwidth))
        
        # include circle that is 1" in radius
        circle_1a = pg.ScatterPlotItem(pxMode=False,size=ang_pix,
                                       pen=pg.mkPen('w'),
                                       brush=pg.mkBrush(255,255, 255, 0))
        spots = [{'pos': [float(margin_here)+1,
                          float(margin_here)+1]}]
        circle_1a.addPoints(spots)
        self.view_img1.addItem(circle_1a)

        img_2 = pg.ImageItem()
        img_2.setImage(np.rot90(image_data2)[::-1])
        img_2.setLookupTable(cubehelix(ch_g,ch_s,ch_r,ch_h))
        img_2.setLevels([level_min,level_max])
        self.view_img2 = self.img2.addViewBox()
        self.view_img2.addItem(img_2)
        self.view_img2.setRange(QtCore.QRectF(0,0,pixwidth,pixwidth))

        level_min_res, level_max_res = vminvmax(image_data3)

        img_3 = pg.ImageItem()
        img_3.setImage(np.rot90(image_data3)[::-1])
        img_3.setLookupTable(cubehelix(ch_g,ch_s,ch_r,ch_h))
        #img_3.setLevels([level_min,level_max])
        img_3.setLevels([level_min_res, level_max_res])
        self.view_img3 = self.img3.addViewBox()
        self.view_img3.addItem(img_3)
        self.view_img3.setRange(QtCore.QRectF(0,0,pixwidth,pixwidth))

        def mouseMovedImg(event):
            pos_add = event.pos()
            margin_here = margin_all_objs[self.list_j[self.j]]
            add_x_pos = pos_add[0]-float(margin_here)+x_pos[self.list_j[self.j]]
            add_y_pos = pos_add[1]-float(margin_here)+y_pos[self.list_j[self.j]]
            self.le_addx.setText(str(int(add_x_pos)))
            self.le_addy.setText(str(int(add_y_pos)))
        img_3.mouseClickEvent = mouseMovedImg

        self.lastClicked = []

        def clicked(plot, points):
           
            h = self.list_j[self.j] # for main object
            self.write_params()
                        
            # untoggle fix params
            for param in xrange(len(params_all_lines.keys())):
                param_here = params_all_lines.keys()[param] 
                if self.cb_params[param_here].isChecked():
                    self.cb_params[param_here].toggle()
    
            for po2 in self.lastClicked:
                po2.resetPen()

            for po in points:
                c = po.data() # for other objects
                self.c = c
                
                po.setPen('b',width=3)
                if po.pos()[0]==float(margin_here)+1: # for main object
                    self.c = 999 # to make sure not to change the main obj
                    self.posx.setText(str(np.round(x_pos[h],2)))
                    self.posy.setText(str(np.round(y_pos[h],2)))
                    for param in xrange(len(params_all_lines.keys())):
                        param_here = params_all_lines.keys()[param] 
                        new_text = params_all_lines[param_here][h]
                        self.param_dict[param_here].setText(str(new_text))

                        # reset fix ceck boxes
                        this_param_fix=fix_other_params[param_here][h]
                        if self.cb_params[param_here].isChecked() and \
                           this_param_fix==1:
                            self.cb_params[param_here].toggle()
                        elif this_param_fix==0:
                            self.cb_params[param_here].toggle()

                    self.cb_fit_each.blockSignals(True)
                    self.cb_inc_each.blockSignals(True)
                    self.cb_fit_each.setChecked(True)
                    self.cb_inc_each.setChecked(True)
                    self.cb_fit_each.blockSignals(False) 
                    self.cb_inc_each.blockSignals(False) 
                                        
                else: # display positions of other objects when clicked

                    for param in xrange(len(params_all_lines.keys())):
                        param_here = params_all_lines.keys()[param] 
                        new_text = params_close_here[h][param_here][c]
                        self.param_dict[param_here].setText(str(new_text))

                        # reset fix ceck boxes
                        this_p_fix=fix_other_params_close_here[h][param_here][c]
                        if self.cb_params[param_here].isChecked() and \
                           this_p_fix==1:
                            self.cb_params[param_here].toggle()
                        elif this_p_fix==0:
                            self.cb_params[param_here].toggle()

                    closex = close_here_x[h][c]
                    closey = close_here_y[h][c]
                    self.posx.setText(str(np.round(closex,2)))
                    self.posy.setText(str(np.round(closey,2)))
                    
                    self.cb_fit_each.blockSignals(True)
                    self.cb_inc_each.blockSignals(True)
                    if fit_other_objects[h][self.c] == False:
                        self.cb_fit_each.setChecked(False)
                    if fit_other_objects[h][self.c] == True:
                        self.cb_fit_each.setChecked(True)

                    if inc_other_objects[h][self.c] == False:
                        self.cb_inc_each.setChecked(False)
                    if inc_other_objects[h][self.c] == True:
                        self.cb_inc_each.setChecked(True) 
                    self.cb_fit_each.blockSignals(False) 
                    self.cb_inc_each.blockSignals(False)                   
                    
            self.lastClicked = points

            self.write_fail_all()
            self.write_red_conf_lam_label()
            if use_for != 'psf':
                self.write_cpart_label()
                self.write_comm_label()

        # mark positions of other objects and make them clickable
        h = self.list_j[self.j]
        for c in xrange(len(close_here_x[h])):
            closex = close_here_x[h][c]
            closey = close_here_y[h][c]
            # brush is completely transparent
            s1 = pg.ScatterPlotItem(size=20, pen=pg.mkPen('r'), 
                                    brush=pg.mkBrush(255, 255, 255, 0))
            spots = [{'pos': [float(margin_here)-float(x_pos[h])+closex+1,
                              float(margin_here)-float(y_pos[h])+closey+1],
                      'data': c}]
            s1.addPoints(spots)
            self.view_img1.addItem(s1)
            ## Make all plots clickable
            s1.sigClicked.connect(clicked)

        # now mark position of object in LSDCat
        s1 = pg.ScatterPlotItem(size=20, pen=pg.mkPen('g'), 
                                brush=pg.mkBrush(255, 255, 255, 0))
        spots = [{'pos': [float(margin_here)+1,
                          float(margin_here)+1],
                  'data': 999}]
        s1.addPoints(spots)
        self.view_img1.addItem(s1)
        s1.sigClicked.connect(clicked)

        # now mark position of other objects in LSDCat
        # This is not really helpful, though...
        #for disp_other_obj in xrange(len(x_pos_all)):
        #    s1 = pg.ScatterPlotItem(size=20, pen=pg.mkPen('g'), 
        #                            brush=pg.mkBrush(255, 255, 255, 0))
        #    spots=[{'pos': [margin_here-x_pos[h]+x_pos_all[disp_other_obj]+1,
        #                    margin_here-y_pos[h]+y_pos_all[disp_other_obj]+1], 
        #              'data': 999}]
        #    s1.setSymbol('d', update=True, dataSet=None, mask=None)
        #    s1.addPoints(spots)
        #    self.view_img1.addItem(s1)

    ############################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
def main():
    app = QtGui.QApplication(sys.argv)
    main = main_Window()
    main.show()    
    app.exec_()

if __name__ == '__main__':
    main()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### PART 4 - read results and create useful output file/plots

def write_useful_output():

    all_lines = []
    all_lines_errors = []
    all_lines_others = [] # for other objects
    all_lines_others_errors = []
    all_sky = []
    bombed = []

    for obj in xrange(len(x_pos)):

        try:
            log_open = open(out_log_dir+'fit_'+band+'_'+str(ID_obj[obj])+'.log','r')
            lines = log_open.readlines()
            # since there are two * in the log file if the fit went wrong
            first = True 
            other_lines = []
            other_lines_errors = []
            for line in xrange(len(lines)):
                if obj_type in lines[line] and first == False:
                    other_lines.append(lines[line])
                    other_lines_errors.append(lines[line+1])
                if obj_type in lines[line] and first == True:
                    all_lines.append(lines[line])
                    all_lines_errors.append(lines[line+1])
                    first = False # the first line is for the main object
                if 'sky' in lines[line]:
                    if 'skycell' in lines[line]:
                        continue
                    all_sky.append(lines[line])
            all_lines_others.append(other_lines)
            all_lines_others_errors.append(other_lines_errors)
            log_open.close()
        except IOError:
            print('   No file for object '+str(ID_obj[obj]))
            bombed_bool[obj] = True
            all_lines.append('bombed')
            all_lines_others.append(['bombed'])
            all_sky.append(['bombed'])


    # read the last parameters of the lines where it worked

    failed = []
    pos_final_x = []
    pos_final_y = []
    pos_final_errors_x = []
    pos_final_errors_y = []

    pos_final_others_x = []
    pos_final_others_y = []
    pos_final_others_error_x = []
    pos_final_others_error_y = []

    pos_final_others_ra = []
    pos_final_others_dec = []

    sky_final = []

    for p in params_all_lines.keys():
        position_here = numbers_for_file[p][1]
        p_array = []
        
        # for main object
        for obj in xrange(len(ID_obj)):
            if ID_obj[obj] not in bombed:
                line = all_lines[obj].split()
                line_errors = all_lines_errors[obj].split()
                value_here = line[position_here]
                value_here_error = line_errors[position_here]

                params_all_lines[p][list_j[obj]] = float(replace(value_here))
                params_all_lines_errors[p][list_j[obj]] = float(replace(value_here_error))

                if '*' in value_here:
                    failed_bool[obj] = True
                else:
                    failed_bool[obj] = False

    # if there were multiple objects, save info for all objects

    params_all_lines_others = copy.deepcopy(params_all_lines)
    params_all_lines_others_errors = copy.deepcopy(params_all_lines_errors)
    
    for p in params_all_lines.keys():
        position_here = numbers_for_file[p][1]
        p_array = []
        params_all_lines_others[p] = [[] for i in list_j]
        params_all_lines_others_errors[p] = [[] for i in list_j]

        # for other objects
        for obj in xrange(len(ID_obj)):
            if ID_obj[obj] not in bombed and match_clear[obj] == 'multi':
                values_here = []
                values_here_error = []
                for other_obj in xrange(len(all_lines_others[obj])):
                    line = all_lines_others[obj][other_obj].split()
                    line_errors=all_lines_others_errors[obj][other_obj].split()
                    value_here = line[position_here]
                    value_here_error = line_errors[position_here]
                    values_here.append(float(replace(value_here)))
                    values_here_error.append(float(replace(value_here_error)))

                params_all_lines_others[p][list_j[obj]]=values_here
                params_all_lines_others_errors[p][list_j[obj]]=values_here_error

    # write positions in array
    for obj in xrange(len(ID_obj)):
        if ID_obj[obj] not in bombed:
            line = all_lines[obj].split()
            line_errors = all_lines_errors[obj].split()

            # if the positions are fixed, they are at a different position
            # in 'line'
            if line[2] !='(':
                pos_final_x.append(float(replace(line[2])))
                pos_final_y.append(float(replace(line[3])))
                pos_final_errors_x.append(float(replace(line_errors[2])))
                pos_final_errors_y.append(float(replace(line_errors[3])))
            else:
                pos_final_x.append(float(replace(line[3])))
                pos_final_y.append(float(replace(line[4])))
                pos_final_errors_x.append(float(replace(line_errors[3])))
                pos_final_errors_y.append(float(replace(line_errors[4])))

    # write positions in array for other objects
    for obj in xrange(len(ID_obj)):
        if ID_obj[obj] not in bombed and match_clear[obj] == 'multi' and \
           len(all_lines_others[obj])>0:
            pos_here_x = []
            pos_here_y = []
            pos_here_error_x = []
            pos_here_error_y = []
            for other_obj in xrange(len(all_lines_others[obj])):
                line = all_lines_others[obj][other_obj].split()
                line_errors = all_lines_others_errors[obj][other_obj].split()
            
                if line[2] !='(': 
                    pos_here_x.append(float(replace(line[2])))
                    pos_here_y.append(float(replace(line[3])))
                    pos_here_error_x.append(float(replace(line_errors[2])))
                    pos_here_error_y.append(float(replace(line_errors[3])))
                else:
                    pos_here_x.append(float(replace(line[3])))
                    pos_here_y.append(float(replace(line[4])))
                    pos_here_error_x.append(float(replace(line_errors[3])))
                    pos_here_error_y.append(float(replace(line_errors[4])))

            # convert x,y positions back to ra, dec
            wcs_obj = wcs.WCS(header_hst, relax=False) 
            coords = wcs_obj.wcs_pix2world(pos_here_x,pos_here_y,0)
            pos_here_others_ra  = coords[0]
            pos_here_others_dec = coords[1]
            
            pos_final_others_ra.append(pos_here_others_ra)
            pos_final_others_dec.append(pos_here_others_dec)
            pos_final_others_x.append(pos_here_x)
            pos_final_others_y.append(pos_here_y)
            pos_final_others_error_x.append(pos_here_error_x)
            pos_final_others_error_y.append(pos_here_error_y)
        else:
            pos_final_others_x.append([])
            pos_final_others_y.append([])
            pos_final_others_error_x.append([])
            pos_final_others_error_y.append([])
            pos_final_others_ra.append([])
            pos_final_others_dec.append([])

    # convert x,y positions back to ra, dec
    wcs_obj = wcs.WCS(header_hst, relax=False) 
    coords = wcs_obj.wcs_pix2world(pos_final_x,pos_final_y,0)
    pos_final_ra  = coords[0]
    pos_final_dec = coords[1]

    fail = []  # make array with info if failed or not
    for obj in xrange(len(ID_obj)):
        if failed_bool[obj] == True:
            fail.append('failed')
        elif bombed_bool[obj] == True:
            fail.append('bombed')
        elif ID_obj[obj] in failed:
            fail.append('looked ok')
        else:
            fail.append('ok')

    # if one parameter of the other objects failed, write in list
    fail_others = []
    for obj in xrange(len(ID_obj)):
        fail_each = []
        if ID_obj[obj] not in bombed and match_clear[obj] == 'multi':
            for a in xrange(len(fit_other_objects[obj])):
                failed = False  
                if fit_other_objects[obj][a] == True:
                    for p in params_all_lines.keys():
                        failed_here = failed_because_close_here[obj][p][a]
                        if failed_here == 1:
                            failed = True
                            fail_each.append(True)  
                    if failed == False:
                        fail_each.append(False)
        fail_others.append(fail_each)

    # Now insert additional objects in lists

    ids = []
    params = []
    params_errors = []
    useful = []
    ra = []
    dec = []
    posx = []
    posy = []
    posx_error = []
    posy_error = []
    guo = []
    match = []
    clump = []
    comm = []
    params_names = []
    for p in params_all_lines.keys():
        params_names.append(p)
        par = []
        par_error = []
        for obj in xrange(len(ID_obj)):
            par.append(params_all_lines[p][obj])
            par_error.append(params_all_lines_errors[p][obj])
            if match_clear[obj] == 'multi':
                list_others = params_all_lines_others[p][obj]
                list_others_error = params_all_lines_others_errors[p][obj]
                for other in xrange(len(list_others)):
                    if inc_other_objects[obj][other] == True:
                        other_value = list_others[other]
                        par.append(other_value)
                        other_value_error = list_others_error[other]
                        par_error.append(other_value_error)

        params.append(par)
        params_errors.append(par_error)


    for obj in xrange(len(ID_obj)):
        ids.append(ID_obj[obj])
        use = fail[obj]
        useful.append(use)
        ra.append(pos_final_ra[obj])
        dec.append(pos_final_dec[obj])
        posx.append(pos_final_x[obj])
        posy.append(pos_final_y[obj])
        posx_error.append(pos_final_errors_x[obj])
        posy_error.append(pos_final_errors_y[obj])
        guo.append(close_here_ID[obj])
        match.append(match_clear[obj])
        if match_clumpy[obj]:
            clump.append(1)
        else:
            clump.append(0)
        comm.append(own_comments[obj])
        sky_final.append(all_sky[obj].split()[4])
        num_i = 0
        if match_clear[obj] == 'multi':
            list_others = params_all_lines_others[p][obj]
            for other in xrange(len(list_others)):
                if inc_other_objects[obj][other] == True:
                    ids.append(ID_obj[obj]+"_"+str(num_i))
                    num_i+=1
                    if fail_others[obj][other] == True:
                        useful.append('failed')
                    else:
                        useful.append('ok')
                    ra.append(pos_final_others_ra[obj][other])
                    dec.append(pos_final_others_dec[obj][other])
                    posx.append(pos_final_others_x[obj][other])
                    posy.append(pos_final_others_y[obj][other])
                    posx_error.append(pos_final_others_error_x[obj][other])
                    posy_error.append(pos_final_others_error_y[obj][other])
                    guo.append(close_here_ID[obj])
                    match.append('-')
                    comm.append('-')
                    clump.append(0)
                    sky_final.append('-')
    
    # Now write output fits file

    col1 = [fits.Column(name='UNIQUE_ID',format='20A',array=ids)]
    cols_params = []
    for p in xrange(len(params_all_lines.keys())):
        col = fits.Column(name=params_names[p].replace(" ","_"),format='E',
                          array=params[p])
        cols_params.append(col)
    cols_params_error = []
    for p in xrange(len(params_all_lines_errors.keys())):
        col = fits.Column(name=params_names[p].replace(" ","_")+"_error",
                          format='E',array=params_errors[p])
        cols_params_error.append(col)
    cols = col1+cols_params+cols_params_error
    cols.append(fits.Column(name='USEFUL',format='20A',array=useful))
    cols.append(fits.Column(name='RA',format='E',array=ra))
    cols.append(fits.Column(name='DEC',format='E',array=dec))
    cols.append(fits.Column(name='GUO_ID',format='E',array=guo))
    cols.append(fits.Column(name='X_pix',format='E',array=posx))
    cols.append(fits.Column(name='Y_pix',format='E',array=posy))
    cols.append(fits.Column(name='X_pix_error',format='E',array=posx_error))
    cols.append(fits.Column(name='Y_pix_error',format='E',array=posy_error))
    cols.append(fits.Column(name='MATCH',format='20A',array=match))
    cols.append(fits.Column(name='CLUMPY',format='B',array=clump))
    cols.append(fits.Column(name='Comment',format='250A',array=comm))
    cols.append(fits.Column(name='Sky',format='50A',array=sky_final))

    tbhdu = fits.BinTableHDU.from_columns(cols)

    current_file_name = inspect.getfile(inspect.currentframe()) 
    prihdr = fits.Header()
    prihdr['SCRIPT'] = current_file_name
    prihdr['DATE'] = str(time.strftime("%d/%m/%Y"))
    prihdr['TIME'] = str(time.strftime("%H:%M:%S"))
    prihdr['IN_IMAGE'] = input_image
    prihdr['IN_CAT'] = input_catalogue
    prihdr['Guo_CAT'] = Guo_catalogue
    prihdr['SIGMA'] = str(sigma_image)
    prihdr['PSF_IMG'] = str(psf_image)
    prihdr['PSF_SAMP'] = str(psf_sampling)
    prihdr['BAD_PIX'] = str(bad_pixel)
    prihdr['CONSTR'] = str(param_constr)
    prihdr['MARGIN'] = str(margin_all_objs[obj])
    prihdr['CONVBOX'] = str(xconv)+','+str(yconv)
    prihdr['FUNC'] = obj_type
    prihdu = fits.PrimaryHDU(header=prihdr)
    thdulist = fits.HDUList([prihdu,tbhdu])
    if save_stuff == True:
        thdulist.writeto(results_final+output_catalogue,overwrite=True)

write_useful_output()

# move output images to different folder
for files in os.listdir(image_dir):
    if 'imgblock' in files:
        shutil.copy2(image_dir+files,result_folder)
        os.remove(image_dir+files)

# move all final results to final result folder
if save_stuff == True and use_for !='psf':
    for i in ID_obj:
        for files in os.listdir(out_log_dir):
            if '.log' in files and i in files:
                # test if dir exists
                if not os.path.isdir(results_final+'logs/'):
                    os.makedirs(results_final+'logs/')
                shutil.copy2(out_log_dir+files,results_final+'logs/')           
        for files in os.listdir(input_params_dir):
            if 'params' in files and i in files:
                #test if dir exists
                if not os.path.isdir(results_final+'input/'):
                    os.makedirs(results_final+'input/')
                shutil.copy2(input_params_dir+files,results_final+'input/')
        for files in os.listdir(result_folder):
            if 'imgblock' in files and i in files:
                #test if dir exists
                if not os.path.isdir(results_final+'imgblocks/'):
                    os.makedirs(results_final+'imgblocks/')
                shutil.copy2(result_folder+files,results_final+'imgblocks/')
    print('   Saved everything!')
else:
    print("   Didn't save results, were you just testing stuff?")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print('\n\n >>> Done! :) <<<')
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
