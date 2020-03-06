
import os,sys, importlib, glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl


sys.path.append('/Users/peterpfleiderer/Projects/regioClim_2020/shaped_analysis')
import shaped_analysis; importlib.reload(shaped_analysis)



COU = shaped_analysis.country_analysis(iso='GHA', working_directory='/home/pepflei/regioClim_2020/cou_data/GHA')
COU.download_shapefile()
COU.load_shapefile()
COU.create_masks(input_file='/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/pr/mon_pr_ECEARTH_rcp45_20060101-21001231_BC.nc4', mask_style='pop2015',add_mask_file = '/home/pepflei/CA/masks/population/population_1990_incrLat.nc', add_mask_name='pop2015', add_mask_var='mask')
COU.create_masks(input_file='/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/pr/mon_pr_ECEARTH_rcp45_20060101-21001231_BC.nc4', mask_style='ppp2005',add_mask_file = '/p/projects/tumble/carls/shared_folder/masks/additional_masks/PPP2005.nc', add_mask_name='ppp2005', add_mask_var='PPP2005')
COU.create_masks(input_file='/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/pr/mon_pr_ECEARTH_rcp45_20060101-21001231_BC.nc4', mask_style='latWeight', lat_weighted=True)

COU.load_mask()

for file_name in glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/pr/mon_pr_*_*_*.nc4'):
    model = file_name.split('_')[-4]
    scenario = file_name.split('_')[-3]
    if scenario == 'historical':
        time_cutoff = ['1950-01-01','2006-01-01']
    else:
        time_cutoff = ['2006-01-01','2100-01-01']
    COU.zoom_data(input_file=file_name,var_name='pr',given_var_name='pr',tags={'scenario':scenario,'experiment':'CORDEX','model':model}, time_cutoff=time_cutoff)

for file_name in glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/tas/mon_tas_*_*_*.nc4'):
    model = file_name.split('_')[-4]
    scenario = file_name.split('_')[-3]
    if scenario == 'historical':
        time_cutoff = ['1950-01-01','2006-01-01']
    else:
        time_cutoff = ['2006-01-01','2100-01-01']
    COU.zoom_data(input_file=file_name,var_name='tas',given_var_name='tas',tags={'scenario':scenario,'experiment':'CORDEX','model':model}, time_cutoff=time_cutoff)

for file_name in glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/TXx/mon_TXx_*_*_*.nc4'):
    model = file_name.split('_')[-4]
    scenario = file_name.split('_')[-3]
    if scenario == 'historical':
        time_cutoff = ['1950-01-01','2006-01-01']
    else:
        time_cutoff = ['2006-01-01','2100-01-01']
    COU.zoom_data(input_file=file_name,var_name='tasmax',given_var_name='TXx',tags={'scenario':scenario,'experiment':'CORDEX','model':model}, time_cutoff=time_cutoff)

for file_name in glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMoutput/monthly/RX1/mon_RX1_*_*_*.nc4'):
    model = file_name.split('_')[-4]
    scenario = file_name.split('_')[-3]
    if scenario == 'historical':
        time_cutoff = ['1950-01-01','2006-01-01']
    else:
        time_cutoff = ['2006-01-01','2100-01-01']
    COU.zoom_data(input_file=file_name,var_name='pr',given_var_name='RX1',tags={'scenario':scenario,'experiment':'CORDEX','model':model}, time_cutoff=time_cutoff)

COU.zoom_data(input_file='/p/projects/tumble/carls/shared_folder/data/EWEMBI/mon_year/mon_tas_EWEMBI_1979-2014.nc4', var_name='tas', tags={'experiment':'EWEMBI','scenario':'historical','model':'EWEMBI'}, time_cutoff=['1950-01-01','2006-01-01'])

COU.zoom_data(input_file='/p/projects/tumble/carls/shared_folder/data/EWEMBI/mon_year/mon_pr_EWEMBI_1979-2014.nc4', var_name='pr', tags={'experiment':'EWEMBI','scenario':'historical','model':'EWEMBI'}, time_cutoff=['1950-01-01','2006-01-01'])

COU.zoom_data(input_file='/p/projects/tumble/carls/shared_folder/data/EWEMBI/mon_year/mon_TXx_EWEMBI_1979-2014.nc4', var_name='tasmax', given_var_name='TXx', tags={'experiment':'EWEMBI','scenario':'historical','model':'EWEMBI'}, time_cutoff=['1950-01-01','2006-01-01'])

COU.zoom_data(input_file='/p/projects/tumble/carls/shared_folder/data/EWEMBI/mon_year/mon_RX1_EWEMBI_1979-2014.nc4', var_name='pr', given_var_name='RX1', tags={'experiment':'EWEMBI','scenario':'historical','model':'EWEMBI'}, time_cutoff=['1950-01-01','2006-01-01'])










#
