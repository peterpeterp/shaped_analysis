
import os,sys, importlib
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl


sys.path.append('/Users/peterpfleiderer/Projects/regioClim_2020/shaped_analysis')
import shaped_analysis; importlib.reload(shaped_analysis)

COU = shaped_analysis.country_analysis(iso='BEN', working_directory='/Users/peterpfleiderer/Projects/regioClim_2020/cou_data/BEN')
COU.load_shapefile()
# COU.identify_grid('/Users/peterpfleiderer/Projects/data/SST/COBE_sst_mon.nc')
# COU.create_masks(input_file='/Users/peterpfleiderer/Projects/data/SST/COBE_sst_mon.nc', var_name='SST', mask_style='pop2015',add_mask_file = '/Users/peterpfleiderer/Projects/data/data_universal/population_2015_incrLat.nc', add_mask_name='pop2015', add_mask_var='mask')
# COU.create_masks(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_vws.nc', mask_style='pop2015',add_mask_file = '/Users/peterpfleiderer/Projects/data/data_universal/population_2015_incrLat.nc', add_mask_name='pop2015', add_mask_var='mask')
COU.load_mask()
# # COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/SST/COBE_sst_mon.nc',var_name='sst',given_var_name='SST',tag='test')
# COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_002_prmsl.nc',var_name='var2',given_var_name='mslp',tags={'scenario':'hist','experiment':'CORDEX','model':'mpi'})
# COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_002_prmsl.nc',var_name='var2',given_var_name='mslp',tags={'scenario':'rcp45','experiment':'CORDEX','model':'mpi'})
# COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_002_prmsl.nc',var_name='var2',given_var_name='mslp',tags={'scenario':'hist','experiment':'CORDEX','model':'had'})
# COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_002_prmsl.nc',var_name='var2',given_var_name='mslp',tags={'scenario':'rcp45','experiment':'CORDEX','model':'had'})


# COU.harmonize_data()
COU.load_data()
# COU.area_average()
COU.load_area_average()



ref = COU.select_data_gridded(time_format='monthly',tags={'scenario':'historical','experiment':'CORDEX','var_name':'pr'})
ref_P = [1980,2000]
ref = ref.loc[:,np.datetime64(str(ref_P[0])+'-01-15'):np.datetime64(str(ref_P[1])+'-01-01')]
monthly = ref.groupby('time.month').mean('time')
ref_seasonal = monthly.loc[:,[2,3,5]].mean(axis=1)

proj = COU.select_data_gridded(time_format='monthly',tags={'scenario':'rcp45','experiment':'CORDEX','var_name':'pr'})
proj_P = [2040,2060]
proj = proj.loc[:,np.datetime64(str(proj_P[0])+'-01-15'):np.datetime64(str(proj_P[1])+'-01-01')]
monthly = proj.groupby('time.month').mean('time')
proj_seasonal = monthly.loc[:,[2,3,5]].mean(axis=1)

diff = proj_seasonal - ref_seasonal
mean_diff = diff.mean('model')

agree = mean_diff.copy()*0
for model in diff.model.values:
    agree += np.sign(mean_diff) == np.sign(diff.loc[model])

agree.values[agree.values>2] = np.nan
agree.values[agree.values<3] = 0.5

plt.close('all')
fig,ax = plt.subplots(nrows=1, figsize=(4,4*(agree.shape[0]/agree.shape[1])**0.5),subplot_kw={'projection': ccrs.PlateCarree()})
ax.coastlines(resolution='10m');
x,y=agree.lon.copy(),agree.lat.copy()
x-=np.diff(x,1)[0]/2.
y-=np.diff(y,1)[0]/2.
x=np.append(x,[x[-1]+np.diff(x,1)[0]])
y=np.append(y,[y[-1]+np.diff(y,1)[0]])
x,y=np.meshgrid(x,y)

im = ax.pcolormesh(x,y,mean_diff, cmap='RdBu_r', transform=ccrs.PlateCarree())
ax.pcolormesh(x,y,agree, cmap='Greys', vmin=0, vmax=1, transform=ccrs.PlateCarree())
for shape in COU._region_polygons.values():
    ax.add_geometries(shape, ccrs.PlateCarree(), edgecolor='k',alpha=1,facecolor='none',linewidth=0.5,zorder=50)
cb = plt.colorbar(im, ax=ax)
cb.set_label('color_label', rotation=90)
plt.savefig(COU._working_dir+'plots/map.png', bbox_inches='tight')






plt.close()
for scenario in ['historical','rcp45']:
    tmp = COU.select_data('monthly',tags={'scenario':scenario,'experiment':'CORDEX','var_name':'tas','region':'BEN','mask_style':'latWeight'})
    tmp = tmp[:,np.all(np.isnan(tmp), axis=0)==False]

    relMon = np.zeros(tmp.time.shape[0], np.bool)
    for mon in [1,3,5]:
        relMon[tmp.time.dt.month.values == mon] = True

    tmp = tmp[:,relMon]

    plt.plot(tmp.year.values,tmp.mean(axis=0), label=scenario)
    plt.fill_between(tmp.year.values,tmp.min(axis=0).values,tmp.max(axis=0).values, alpha=0.2)
plt.legend()
plt.ylabel('temperature [K]')
plt.xlabel('year')
plt.savefig(COU._working_dir+'plots/transient.png')









plt.close()
for mask_style,color in zip(['latWeight','pop2015','ppp2005'],['green','orange','red']):
    for scenario in ['historical','rcp45']:
        tmp = COU.select_data('monthly',tags={'scenario':scenario,'experiment':'CORDEX','var_name':'tas','region':'BEN','mask_style':mask_style})
        tmp = tmp[:,np.all(np.isnan(tmp), axis=0)==False]
        tmp = tmp.groupby('time.year').mean('time')
        plt.plot(tmp.year.values,tmp.mean(axis=0), color = color, label=mask_style)
plt.legend()
plt.ylabel('temperature [K]')
plt.xlabel('year')
plt.savefig(COU._working_dir+'plots/different_masks.png')



asdas

# COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_002_prmsl.nc',var_name='var2',given_var_name='mslp',tag='rcp85')
# COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_vws.nc',var_name='var33',given_var_name='VWS',tag='JRA55')






# COU = country_analysis(iso='BEN', working_directory='/home/pepflei/regioClim_2020/cou_data/BEN')
# # COU.download_shapefile()
# COU.load_shapefile()
# # COU.create_masks(input_file='/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_ECEARTH_rcp45_.nc4', mask_style='pop2015',add_mask_file = '/home/pepflei/CA/masks/population/population_1990_incrLat.nc', add_mask_name='pop2015', add_mask_var='mask')
# # COU.create_masks(input_file='/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_ECEARTH_rcp45_.nc4', mask_style='ppp2005',add_mask_file = '/p/projects/tumble/carls/shared_folder/masks/additional_masks/PPP2005.nc', add_mask_name='ppp2005', add_mask_var='PPP2005')
# # COU.create_masks(input_file='/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_ECEARTH_rcp45_.nc4', mask_style='latWeight', lat_weighted=True)
# COU.load_mask()
#
# for file_name in glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_*_*_.nc4'):
#     model = file_name.split('_')[-3]
#     scenario = file_name.split('_')[-2]
#     if scenario == 'historical':
#         time_cutoff = [np.datetime64('1950-01-01'),np.datetime64('2006-01-01')]
#     else:
#         time_cutoff = [np.datetime64('2006-01-01'),np.datetime64('2100-01-01')]
#     COU.zoom_data(input_file=file_name,var_name='pr',given_var_name='pr',tags={'scenario':scenario,'experiment':'CORDEX','model':model}, time_cutoff=time_cutoff)
#
# for file_name in glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/tas/mon_tas_*_*_.nc4'):
#     model = file_name.split('_')[-3]
#     scenario = file_name.split('_')[-2]
#     if scenario == 'historical':
#         time_cutoff = [np.datetime64('1950-01-01'),np.datetime64('2006-01-01')]
#     else:
#         time_cutoff = [np.datetime64('2006-01-01'),np.datetime64('2100-01-01')]
#     COU.zoom_data(input_file=file_name,var_name='tas',given_var_name='tas',tags={'scenario':scenario,'experiment':'CORDEX','model':model}, time_cutoff=time_cutoff)
#
#
# COU.load_data()
# COU.area_average()
# asdas