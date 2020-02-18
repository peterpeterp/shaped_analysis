# -*- coding: utf-8 -*-
'''
Class to analyze climate data on national (& sub national) scale
'''
# Author: Peter Pfleiderer <peter.pfleiderer@climateanalytics.org>
#
# License: GNU General Public License v3.0

import sys,glob,os
import numpy as np
import xarray as xr
import pandas as pd
from shapely.geometry import mapping, Polygon, MultiPolygon, asShape
from shapely.ops import cascaded_union, unary_union
import matplotlib.pylab as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shapereader

class country_analysis(object):

    def __init__(self, iso, working_directory):
        '''
        Prepare directories and meta-data
        iso: str: Country name or iso
        working_directory: path: Path where files will be stored
        seasons: dict: seasons relevant for the country. 'season name':{months in season as int 1-12}
        additional_tag: str: when specified raw data will be stored in a separate directory with additional tag name.
        '''
        self._iso=iso
        self._working_dir=working_directory+'/'
        self._support_dir= '/'+'/'.join(working_directory.split('/')[:-1])+'/support/'

        if os.path.isdir(self._working_dir)==False: os.system('mkdir '+self._working_dir)
        if os.path.isdir(self._working_dir+'/masks')==False: os.system('mkdir '+self._working_dir+'/masks')
        if os.path.isdir(self._working_dir+'/plots')==False: os.system('mkdir '+self._working_dir+'/plots')
        if os.path.isdir(self._working_dir+'/gridded_data')==False:os.system('mkdir '+self._working_dir+'/gridded_data')
        if os.path.isdir(self._working_dir+'/processed_data')==False:os.system('mkdir '+self._working_dir+'/processed_data')

        self._masks={}
        self._region_names = {}
        self._masks = {}
        self._grid_dict = {}

    def download_shapefile(self):
        '''
        to be tested
        '''
        current_dir=os.getcwd()
        os.chdir(self._working_dir)
        os.system('mkdir '+self._iso+'_adm_shp')
        os.system('wget biogeo.ucdavis.edu/data/gadm2.8/shp/'+self._iso+'_adm_shp.zip')
        os.chdir(self._working_dir + self._iso+'_adm_shp')
        os.system('unzip ../'+self._iso+'_adm_shp.zip')
        os.chdir(current_dir)

    def load_shapefile(self):

        adm_shapefiles=shapereader.Reader(self._working_dir+self._iso+'_adm_shp/'+self._iso+'_adm1.shp').records()
        # collect all shapes of region
        self._region_polygons={}
        for item in adm_shapefiles:
            shape,region=item.geometry,item.attributes
            region = {k.lower():v for k,v in region.items()}
            name_full = u''+region['name_1']
            try:
                name=u''+region[u'hasc_1']
            except:
                print(region)
                name=u''+region[u'type_1']
            self._region_names[name]=name_full
            # simplify could be added here to speed up things
            try:
                self._region_polygons[name]=MultiPolygon(shape)
            except:
                self._region_polygons[name]=Polygon(shape)

        name=self._iso
        self._region_names[name]=name
        try:
            adm_shapefiles=shapereader.Reader(self._working_dir+self._iso+'_adm_shp/'+self._iso+'_adm0.shp').records()
            self._region_polygons[self._iso]=MultiPolygon(next(adm_shapefiles).geometry)
        except:
            adm_shapefiles=shapereader.Reader(self._working_dir+self._iso+'_adm_shp/'+self._iso+'_adm0.shp').records()
            self._region_polygons[self._iso]=Polygon(next(adm_shapefiles).geometry)

        print(self._region_names)

    ##########
    # MASKS
    ##########

    def load_mask(self):
        for mask_file in glob.glob(self._working_dir+'/masks/*.nc4'):
            small_grid = mask_file.split('_')[-2]
            if small_grid not in self._masks.keys():	self._masks[small_grid]={}

            mask_style = mask_file.split('_')[-1].split('.')[0]
            nc_mask = xr.open_dataset(mask_file)
            self._masks[small_grid][mask_style] = nc_mask['mask']
            self._grid_dict[nc_mask.attrs['original_grid']] = small_grid

    def identify_grid(self,input_file,lat_name='lat',lon_name='lon'):
        '''
        get information about grid of input data
        input_file: file_path: file to be analyzed
        lat_name: str: name of latitude variable
        lon_name: str: name of longitude variable
        '''
        nc_in = xr.open_dataset(input_file)
        lat = nc_in[lat_name][:]
        lon = nc_in[lon_name][:]

        if len(lat.shape)==2:
            lat=lat[:,0]
            lon=lon[0,:]

        if max(lon)>200:	lon_shift=-180.0
        else:				lon_shift=0.0
        lon+=lon_shift

        nx = len(lon)	;	ny = len(lat)
        grid=str(str(ny)+'x'+str(nx)+'lat'+str(lat.values[0])+'to'+str(lat.values[-1])+'lon'+str(lon.values[0])+'to'+str(lon.values[-1])).replace('.','p')
        nc_in.close()

        # self._raw_lon = lon.values
        # self._raw_lat = lat.values
        # self._raw_grid = grid
        # slef._raw_lon_shift = lon_shift

        return lon.values,lat.values,grid,lon_shift

    def regrid_additional_mask(self,add_mask_file,add_mask_name,add_mask_var,grid,lon,lat,shift):
        '''
        regrid population masks
        grid: str: name of the grid
        lon: array: longitudes
        lat: array: latitudes
        shift: int: number of elements to roll around longitude axis. Has to be considered since the longitude axis might be shifted (see get_grid_polygons)
        add_mask_file: file_path: path of used population mask
        add_mask_name: str: name of the created mask
        '''
        regridded_mask_file = self._support_dir+add_mask_name+'_'+grid+'.nc'
        if os.path.isfile(regridded_mask_file) == False:
            mygrid=open(self._support_dir+grid+'.txt','w')
            mygrid.write('gridtype=lonlat\nxsize='+str(len(lon))+'\nysize='+str(len(lat))+'\nxfirst='+str(lon[0])+'\nxinc='+str(np.mean(np.diff(lon,1)))+'\nyfirst='+str(lat[0])+'\nyinc='+str(np.mean(np.diff(lat,1))))
            mygrid.close()
            os.system('cdo remapbil,'+self._support_dir+grid+'.txt '+add_mask_file+' '+regridded_mask_file)

        add_mask = xr.open_dataset(regridded_mask_file)[add_mask_var].squeeze()
        add_mask = np.roll(add_mask.values,shift,axis=1)
        return add_mask

    def get_grid_polygons(self,grid,lon,lat,lon_shift):
        '''
        create polygons for each grid-cell
        grid: str: name of the grid
        lon: array: longitudes
        lat: array: latitudes
        lon_shift: float: deg longitudes that have to be shifted to be on a -180 to 180 grid (computed in identify_grid)
        '''
        # loop over the grid to get grid polygons
        nx = len(lon)	;	ny = len(lat)

        grid_polygons = np.empty((nx,ny),dtype=Polygon)
        dx = np.zeros((nx))
        dy = np.zeros((ny))
        dx[1:] = np.abs(np.diff(lon,1))
        dx[0] = dx[1]
        dy[1:] = np.abs(np.diff(lat,1))
        dy[0] = dy[1]
        for i in range(nx):
            x1 = lon[i]-dx[i]/2.
            x2 = lon[i]+dx[i]/2.
            for j in range(ny):
                y1 = lat[j]-dy[j]/2.
                y2 = lat[j]+dy[j]/2.
                grid_polygons[i,j] = Polygon([(x1,y1),(x1,y2),(x2,y2),(x2,y1)])
                #grid_polygons[i,j] = Polygon([(y1,x1),(y1,x2),(y2,x2),(y2,x1)])

        # since the lon axis has been shifted, masks and outputs will have to be shifted as well. This shift is computed here
        lon-=lon_shift
        shift = len(lon)-np.where(lon==lon[0]-lon_shift)[0][0]

        self._masks[grid]['lat_mask']=lat
        self._masks[grid]['lon_mask']=lon

        return grid_polygons,shift

    def grid_polygon_overlap(self,grid,lon,lat,grid_polygons,country_polygons,shift,mask_style,ext_poly,add_mask=None):
        '''
        Compute overlap betwwen grid polygons (get_grid_polygons) and country polygons
        grid: str: name of the grid
        lon: array: longitudes
        lat: array: latitudes
        grid_polygons: list: List of polygons created in get_grid_polygons
        country_polgons: list: list of polygons representing the country
        shift: int: number of elements to roll around longitude axis. Has to be considered since the longitude axis might be shifted (see get_grid_polygons)
        mask_style: str: Can be 'lat_weighted' or population weighted. If population weighted, mask_style is a given name
        est_poly: Polygon: Polygon limiting the are where overlaps are computed
        name: str: country-name or region name
        add_mask: np.array: population mask from regrid_add_mask
        '''
        nx = len(lon)	;	ny = len(lat)

        overlap = np.zeros((ny,nx))
        for i in range(nx):
            for j in range(ny):
                # check gridcell is relevant
                if grid_polygons[i,j].intersects(ext_poly):
                    # get fraction of grid-cell covered by polygon
                    intersect = grid_polygons[i,j].intersection(country_polygons).area/grid_polygons[i,j].area*country_polygons.area
                    if add_mask is not None:
                        # population weighting
                        overlap[j,i] = intersect*add_mask[j,i]
                    if mask_style=='lat_weighted':
                        # multiply overlap with latitude weighting
                        overlap[j,i] = intersect*np.cos(np.radians(lat[j]))

        # renormalize overlap to get sum(mask)=1
        overlap_sum=sum(overlap.copy().flatten())
        if overlap_sum!=0:
            output=overlap.copy()/overlap_sum
            # mask zeros
            output[output==0]=np.nan
            output=np.ma.masked_invalid(output)
            # shift back to original longitudes
            return np.roll(output,shift,axis=1)
        else:
            print('something went wrong with the mask')
            return False

    def create_masks(self,input_file,mask_style,lat_weighted=True,add_mask_file='',add_mask_name='add_mask',add_mask_var='mask',lat_name='lat',lon_name='lon'):
        '''
        create country mask
        input_file: str: location of example input data (required for the identification of the grid)
        shape_file: str: location of the shape_file used to identify country borders
        mask_style: str: name under which the mask will be stored (important for further analysis)
        add_mask_file: str: location of population mask (netcdf file) used for population weighted country mask
        overwrite: bool: if True, old files is deleted, new mask created
        lat_name: str: name of latitude variable in netcdf file
        lon_name: str: name of longitude variable in netcdf file
        '''

        lon,lat,grid,lon_shift = self.identify_grid(input_file,lat_name,lon_name)

        if grid not in self._masks.keys():
            self._masks[grid]={}

        grid_polygons,shift = self.get_grid_polygons(grid,lon,lat,lon_shift)

        country_polygons = self._region_polygons[self._iso]

        # get boundaries for faster computation
        x1, y1, x2, y2 = country_polygons.bounds
        xmin, xmax, ymin, ymax = min([x1,x2]), max([x1,x2]), min([y1,y2]), max([y1,y2])
        ext = [(xmin,ymin),(xmin,ymax),(xmax,ymax),(xmax,ymin),(xmin,ymin)]
        ext_poly = Polygon(ext)

        # load population mask
        if add_mask_file=='':
            add_mask = np.ones((len(lat),len(lon)))
        else:
            add_mask = self.regrid_additional_mask(add_mask_file,add_mask_name,add_mask_var,grid,lon,lat,shift)

        # preprare output
        out_mask = xr.DataArray(data=np.zeros([len(self._region_polygons.keys()),len(lat),len(lon)])*np.nan, coords={'region':list(self._region_polygons.keys()),'lat':lat,'lon':lon}, dims=['region','lat','lon'])

        # compute overlap
        for region_name,polygons in self._region_polygons.items():
            out_mask.loc[region_name,:,:] = self.grid_polygon_overlap(grid, lon, lat, grid_polygons, polygons, shift, mask_style, ext_poly, add_mask)

        # get relevant extend
        cou_mask=np.ma.getdata(out_mask.loc[self._iso])
        lons=sorted(np.where(np.isfinite(np.nanmean(cou_mask,0)))[0])
        lon_=lon[lons[0]:lons[-1]+1]
        lats=sorted(np.where(np.isfinite(np.nanmean(cou_mask,1)))[0])
        lat_=lat[lats[0]:lats[-1]+1]

        small_grid=str(str(len(lat_))+'x'+str(len(lon_))+'lat'+str(lat_[0])+'to'+str(lat_[-1])+'lon'+str(lon_[0])+'to'+str(lon_[-1])).replace('.','p')
        self._grid_dict[grid] = small_grid
        if small_grid not in self._masks.keys():	self._masks[small_grid]={}

        zoomed_mask = out_mask[:,lats[0]:lats[-1]+1,lons[0]:lons[-1]+1]
        self._masks[small_grid][mask_style] = zoomed_mask
        mask_file=self._working_dir+'/masks/'+self._iso+'_'+small_grid+'_'+mask_style+'.nc4'
        nc_out_mask = xr.Dataset({'mask': zoomed_mask})
        nc_out_mask.attrs['original_grid'] = grid
        nc_out_mask.attrs['grid'] = small_grid
        nc_out_mask.attrs['mask_style'] = mask_style
        nc_out_mask.to_netcdf(self._working_dir+'/masks/'+self._iso+'_'+small_grid+'_'+mask_style+'.nc4')

    def zoom_data(self,input_file,var_name,given_var_name=None,tags={},lat_name='lat',lon_name='lon'):
        lon,lat,grid,lon_shift = self.identify_grid(input_file,lat_name,lon_name)
        input = xr.open_dataset(input_file)[var_name].squeeze()
        mask = self._masks[self._grid_dict[grid]][list(self._masks[self._grid_dict[grid]].keys())[0]]

        zoomed_data = input.loc[:,mask.lat,mask.lon]
        if given_var_name is None:
            given_var_name = var_name
        out_data = xr.Dataset({given_var_name: zoomed_data})
        out_data.attrs['original_grid'] = grid
        out_data.attrs['grid'] = self._grid_dict[grid]
        out_data.attrs['var_name_original'] = var_name
        tag = ''
        for key,val in tags.items():
            out_data[given_var_name].attrs['tag_'+key] = val
            tag += '_'+key+'-'+val
        out_data.to_netcdf(self._working_dir+'/gridded_data/'+self._iso+'_'+self._grid_dict[grid]+'_'+given_var_name+tag+'.nc4')

    def load_data(self):
        # get information about all files
        grids = [nn.split('_')[-2] for nn in glob.glob(self._working_dir+'/masks/*.nc4')]
        tag_dict = {grid:{} for grid in grids}
        coords_dict = {grid:{} for grid in grids}
        for file_name in glob.glob(self._working_dir+'/gridded_data/*.nc4'):
            nc_in = xr.open_dataset(file_name)
            var_name = file_name.split('/')[-1].split('_')[2]
            grid = nc_in.attrs['grid']
            for key,val in nc_in[var_name].attrs.items():
                if 'tag' in key:
                    if key.split('_')[-1] not in tag_dict[grid].keys():
                        tag_dict[grid][key.split('_')[-1]] = []
                    tag_dict[grid][key.split('_')[-1]].append(val)
            if 'time' not in coords_dict[grid].keys():
                coords_dict[grid]['time'] = nc_in['time']
            else:
                coords_dict[grid]['time'].combine_first(nc_in['time'])
            if 'lon' not in coords_dict[grid].keys():
                coords_dict[grid]['lon'] = nc_in['lon']
            if 'lat' not in coords_dict[grid].keys():
                coords_dict[grid]['lat'] = nc_in['lat']

        # create data arrays
        self._data = {}
        for grid,tmp_dict in tag_dict.items():
            for key,val in tmp_dict.items():
                tmp_dict[key] = np.unique(val)
            coords = tmp_dict.copy()
            coords.update({key:val.values for key,val in coords_dict[grid].items()})
            dims = ['time','lat','lon']
            for key in sorted(coords.keys()):
                if key not in dims:
                    dims = [key] + dims

            if len(list(coords.keys())) > 0:
                self._data[grid] = xr.DataArray(data=np.zeros([len(coords[key]) for key in dims])*np.nan, coords=coords, dims=dims)

        # fill data arrays
        for file_name in glob.glob(self._working_dir+'/gridded_data/*.nc4'):
            nc_in = xr.open_dataset(file_name)
            var_name = file_name.split('/')[-1].split('_')[2]
            grid = nc_in.attrs['grid']

            indices = [nc_in[var_name].attrs['tag_'+dim_] for dim_ in self._data[grid].dims[:-3]]

            '''
            this is really ugly but I don't know how to make this differently
            '''
            if len(indices) == 0:
                self._data[grid].loc[nc_in.time,nc_in.lat,nc_in.lon]
            if len(indices) == 1:
                self._data[grid].loc[indices[0],nc_in.time,nc_in.lat,nc_in.lon]
            if len(indices) == 2:
                self._data[grid].loc[indices[0],indices[1],nc_in.time,nc_in.lat,nc_in.lon]
            if len(indices) == 3:
                self._data[grid].loc[indices[0],indices[1],indices[2],nc_in.time,nc_in.lat,nc_in.lon] = nc_in[var_name]
            if len(indices) == 4:
                self._data[grid].loc[indices[0],indices[1],indices[2],indices[3],nc_in.time,nc_in.lat,nc_in.lon] = nc_in[var_name]

        print('Loaded data into self._data')
        for key, val in self._data.items():
            print('on grid: '+key)
            print(val.coords,'\n')







# COU = country_analysis(iso='BEN', working_directory='/home/pepflei/regioClim_2020/cou_data/BEN')
# # COU.download_shapefile()
# # COU.load_shapefile()
# # COU.create_masks(input_file='/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_ECEARTH_rcp45_.nc4', mask_style='pop2015',add_mask_file = '/home/pepflei/CA/masks/population/population_1990_incrLat.nc', add_mask_name='pop2015', add_mask_var='mask')
# COU.load_mask()
# COU.zoom_data(input_file='/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_ECEARTH_rcp45_.nc4',var_name='pr',given_var_name='pr',tag='rcp45')
# COU.zoom_data(input_file='/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_ECEARTH_historical_.nc4',var_name='pr',given_var_name='pr',tag='rcp45')
#
# COU.load_data()
#
# asdas


COU = country_analysis(iso='BEN', working_directory='/Users/peterpfleiderer/Projects/regioClim_2020/country_analysis/data/BEN')
# COU.load_shapefile()
# COU.identify_grid('/Users/peterpfleiderer/Projects/data/SST/COBE_sst_mon.nc')
# COU.regrid_additional_mask(add_mask_file='/Users/peterpfleiderer/Projects/data/data_universal/population_2015_incrLat.nc', mask_name='pop2015', input_file='/Users/peterpfleiderer/Projects/data/SST/COBE_sst_mon.nc')
# COU.create_masks(input_file='/Users/peterpfleiderer/Projects/data/SST/COBE_sst_mon.nc', var_name='SST', mask_style='pop2015',add_mask_file = '/Users/peterpfleiderer/Projects/data/data_universal/population_2015_incrLat.nc', add_mask_name='pop2015', add_mask_var='mask')
# COU.create_masks(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_vws.nc', mask_style='pop2015',add_mask_file = '/Users/peterpfleiderer/Projects/data/data_universal/population_2015_incrLat.nc', add_mask_name='pop2015', add_mask_var='mask')
# COU.load_mask()
# # COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/SST/COBE_sst_mon.nc',var_name='sst',given_var_name='SST',tag='test')
# COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_002_prmsl.nc',var_name='var2',given_var_name='mslp',tags={'scenario':'hist','experiment':'CORDEX','model':'mpi'})
# COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_002_prmsl.nc',var_name='var2',given_var_name='mslp',tags={'scenario':'rcp45','experiment':'CORDEX','model':'mpi'})
# COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_002_prmsl.nc',var_name='var2',given_var_name='mslp',tags={'scenario':'hist','experiment':'CORDEX','model':'had'})
# COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_002_prmsl.nc',var_name='var2',given_var_name='mslp',tags={'scenario':'rcp45','experiment':'CORDEX','model':'had'})


COU.load_data()

# COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_002_prmsl.nc',var_name='var2',given_var_name='mslp',tag='rcp85')
# COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_vws.nc',var_name='var33',given_var_name='VWS',tag='JRA55')



#
