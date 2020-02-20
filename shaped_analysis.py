# -*- coding: utf-8 -*-
'''
Class to analyze climate data on national (& sub national) scale
'''
# Author: Peter Pfleiderer <peter.pfleiderer@climateanalytics.org>
#
# License: GNU General Public License v3.0

import sys,glob,os,pickle
from collections import Counter
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

def save_pkl(obj, name ):
	with open(name, 'wb') as f:
		pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_pkl(name ):
	with open( name, 'rb') as f:
		return pickle.load(f)

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
            self._masks[small_grid][mask_style] = nc_mask['mask'].sortby(['lat', 'lon'])

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
        out_data.to_netcdf(self._working_dir+'/raw_data/'+self._iso+'_'+self._grid_dict[grid]+'_'+given_var_name+tag+'.nc4')

    def harmonize_time(self,tmp_time):
        _, index = np.unique(tmp_time, return_index=True)
        if np.median(np.diff(tmp_time[index].dt.year,1)[1:]) == 1: # yearly
            tmp_time = np.array([np.datetime64(str(dd)[:4]+'-06-15T00:00:00.000000000') for dd in tmp_time.values])
            time_format = 'yearly'
        elif np.median(np.diff(tmp_time[index].dt.month,1)[1:]) == 1: # monthly
            tmp_time = np.array([np.datetime64(str(dd)[:8]+'15T00:00:00.000000000') for dd in tmp_time.values])
            time_format = 'monthly'
        elif np.median(np.diff(tmp_time[index].dt.day,1)[1:]) == 1: # monthly
            tmp_time = np.array([np.datetime64(str(dd)[:10]+'T00:00:00.000000000') for dd in tmp_time.values])
            time_format = 'daily'
        else:
            time_format = None
            print('time format not yearly, not monthly, not daily -> unknown')

        return tmp_time,time_format, index

    def harmonize_data(self):
        # get information about all files
        grids = [nn.split('_')[-2] for nn in glob.glob(self._working_dir+'/masks/*.nc4')]
        self._tags = ['var_name']
        self._tag_combinations = {grid:{} for grid in grids}
        coords_dict = {grid:{'monthly':{'var_name':[]},'yearly':{'var_name':[]},'daily':{'var_name':[]}} for grid in grids}
        coords_dict_area = {'monthly':{'var_name':[]},'yearly':{'var_name':[]},'daily':{'var_name':[]}}
        for file_name in glob.glob(self._working_dir+'/raw_data/*.nc4'):
            nc_in = xr.open_dataset(file_name)
            var_name = file_name.split('/')[-1].split('_')[2]
            grid = nc_in.attrs['grid']

            # identify and harmonize time format
            tmp_time,time_format,index = self.harmonize_time(nc_in['time'])
            if time_format is None:
                print(file_name)
                asdas

            # store coordinates
            if time_format is not None:
                tmp_time = xr.DataArray(tmp_time, coords={'time':tmp_time}, dims=['time'])
                if 'time' not in coords_dict[grid][time_format].keys():
                    coords_dict[grid][time_format]['time'] = tmp_time
                    coords_dict_area[time_format]['time'] = nc_in['time']
                else:
                    coords_dict[grid][time_format]['time'] = xr.concat([coords_dict[grid][time_format]['time'],tmp_time], dim='time')
                    coords_dict_area[time_format]['time'] = xr.concat([coords_dict_area[time_format]['time'],tmp_time], dim='time')
                if 'lon' not in coords_dict[grid][time_format].keys():
                    coords_dict[grid][time_format]['lon'] = nc_in['lon']
                    coords_dict_area[time_format]['lon'] = nc_in['lon']
                if 'lat' not in coords_dict[grid][time_format].keys():
                    coords_dict[grid][time_format]['lat'] = nc_in['lat']
                    coords_dict_area[time_format]['lat'] = nc_in['lat']

                # store tags
                tag_combi = {'var_name':var_name}
                for key,val in nc_in[var_name].attrs.items():
                    if 'tag' in key:
                        tag_combi[key.split('_')[-1]] = val
                        if key.split('_')[-1] not in coords_dict[grid][time_format].keys():
                            coords_dict[grid][time_format][key.split('_')[-1]] = []
                        if key.split('_')[-1] not in self._tags:
                            self._tags += [key.split('_')[-1]]
                        coords_dict[grid][time_format][key.split('_')[-1]].append(val)
                if var_name not in coords_dict[grid][time_format]:
                    coords_dict[grid][time_format]['var_name'] += [var_name]

                for key,val in nc_in[var_name].attrs.items():
                    if 'tag' in key:
                        tag_combi[key.split('_')[-1]] = val
                        if key.split('_')[-1] not in coords_dict_area[time_format].keys():
                            coords_dict_area[time_format][key.split('_')[-1]] = []
                        if key.split('_')[-1] not in self._tags:
                            self._tags += [key.split('_')[-1]]
                        coords_dict_area[time_format][key.split('_')[-1]].append(val)
                if var_name not in coords_dict_area[time_format]:
                    coords_dict_area[time_format]['var_name'] += [var_name]

                if time_format not in self._tag_combinations[grid].keys():
                    self._tag_combinations[grid][time_format] = {}
                self._tag_combinations[grid][time_format][file_name] = tag_combi

        save_pkl(self._tag_combinations, self._working_dir+'meta_info.pkl')
        for time_format,coords in coords_dict_area.items():
            for key,val in coords.items():
                coords[key] = np.unique(val)
        save_pkl(coords_dict_area, self._working_dir+'coords_areaAverage.pkl')

        # create data arrays
        self._data = {grid:{} for grid in grids}
        for grid,tmp_dict in coords_dict.items():
            for time_format,coords in tmp_dict.items():
                if len(coords['var_name']) > 0:
                    for key,val in coords.items():
                        coords[key] = np.unique(val)

            dims = ['time','lat','lon']
            for key in sorted(coords.keys()):
                if key not in dims:
                    dims = [key] + dims

            coords = {dim:coords[dim] for dim in dims}
            if len(dims) > 0:
                self._data[grid][time_format] = xr.DataArray(data=np.zeros([len(coords[key]) for key in dims])*np.nan, coords=coords, dims=dims)

        # fill data arrays
        for grid,tmp1_dict in self._tag_combinations.items():
            for time_format, tmp2_dict in tmp1_dict.items():
                for file_name, tag_dict in tmp2_dict.items():
                    nc_in = xr.open_dataset(file_name)
                    var_name = file_name.split('/')[-1].split('_')[2]
                    tmp_time,time_format,index = self.harmonize_time(nc_in['time'])

                    nc_in['time'] = tmp_time
                    nc_in = nc_in.isel(time=index)

                    indices = [tag_dict[dim] for dim  in self._data[grid][time_format].dims if dim not in ['time','lat','lon']]

                    '''
                    this is really ugly but I don't know how to make this differently
                    '''
                    if len(indices) == 0:
                        self._data[grid][time_format].loc[nc_in.time,nc_in.lat,nc_in.lon]
                    if len(indices) == 1:
                        self._data[grid][time_format].loc[indices[0],nc_in.time,nc_in.lat,nc_in.lon]
                    if len(indices) == 2:
                        self._data[grid][time_format].loc[indices[0],indices[1],nc_in.time,nc_in.lat,nc_in.lon]
                    if len(indices) == 3:
                        self._data[grid][time_format].loc[indices[0],indices[1],indices[2],nc_in.time,nc_in.lat,nc_in.lon] = nc_in[var_name]
                    if len(indices) == 4:
                        self._data[grid][time_format].loc[indices[0],indices[1],indices[2],indices[3],nc_in.time,nc_in.lat,nc_in.lon] = nc_in[var_name]

                if len(list(tmp2_dict.keys())) > 0:
                    self._data[grid][time_format].sortby(['lat', 'lon','time'])

        print('Loaded data into self._data')
        for grid, tmp1_dict in self._data.items():
            for time_format, tmp in tmp1_dict.items():
                print('grid: '+grid)
                print('time format: '+time_format)
                print(tmp.coords,'\n')
                xr.Dataset({'data':tmp}).to_netcdf(self._working_dir+'gridded_data/'+'_'.join([self._iso,grid,time_format])+'.nc')

    def load_data(self):
        self._data = {}
        for file_name in glob.glob(self._working_dir+'/gridded_data/*.nc'):
            nc_in = xr.open_dataset(file_name)
            time_format = file_name.split('_')[-1].split('.')[0]
            grid = file_name.split('_')[-2]
            if grid not in self._data.keys():
                self._data[grid] = {}
            self._data[grid][time_format] = nc_in['data']

    def area_average(self, regions=None):
        self._tag_combinations = load_pkl(self._working_dir+'meta_info.pkl')
        coords_area = load_pkl(self._working_dir+'coords_areaAverage.pkl')
        self._areaAverage = {}
        time_formats = np.array([[time_format for time_format in tmp1.keys()] for tmp1 in self._tag_combinations.values()]).flatten()
        for time_format in time_formats:
            dims = [dim for dim in coords_area[time_format].keys() if dim not in ['lat','lon','time']] + ['time']
            coords = {dim:coords_area[time_format][dim]  for dim in dims}
            self._areaAverage[time_format] = xr.DataArray(data=np.zeros([len(coords[key]) for key in dims])*np.nan, coords=coords, dims=dims)

        for region in self._region_names.keys():
            for time_format in time_formats:
                all_tags = []
                for grid,tmp1 in self._tag_combinations.items():
                    for time_format,tmp2 in tmp1.items():
                        all_tags += list(tmp2[list(tmp2.keys())[0]].keys())

                csv = pd.DataFrame(columns=all_tags+['mask_style','time','value'])

            for grid,tmp1 in self._tag_combinations.items():
                for time_format,tmp2 in tmp1.items():
                    for tag_combi in tmp2.values():
                        tmp = self._data[grid][time_format].loc[tuple(tag_combi[dim] for dim in self._data[grid][time_format].dims[:-3])].copy()

                        for mask_style,mask in self._masks[grid].items():
                            tmp *= mask.loc[region]

                            av = tmp.sum(axis=(-2,-1))

                            dims = self._data[grid][time_format].dims[:-2]
                            coords = {dim:np.array([tag_dict[dim]])  for dim in dims if dim != 'time'}
                            coords['time'] = self._data[grid][time_format].time
                            slice_ = av.copy().values
                            for i in range(len(dims)-1):
                                slice_ = np.expand_dims(slice_,0)
                            data = xr.DataArray(data = slice_, coords=coords, dims=dims)
                            if time_format not in self._areaAverage.keys():
                                self._areaAverage[time_format] = data
                            else:
                                self._areaAverage[time_format] = xr.auto_combine([self._areaAverage[time_format],data])


                            print(self._areaAverage[time_format].shape)

                            available = np.where(np.all(np.isnan(tmp.values), axis=(1,2))==False)[0]
                            av = av[available]
                            tmp_ = xr.Dataset({'value':av}).to_dataframe()
                            tmp_ = tmp_.drop(columns=['region'])
                            tmp_.reset_index(inplace=True)
                            tmp_['mask_style'] = mask_style

                            csv = csv.append(tmp_, sort=True)

            for time_format in time_formats:
                csv.to_csv(self._working_dir+'areaAverage/'+region+'_'+time_format+'_areaAverage.csv')

        # self._areaAverage = {}
        # for time_format in time_formats:
        #     self._areaAverage[time_format] = []
        # for grid,tmp1 in self._data.items():
        #     for time_format,data in tmp1.items():
        #         # if time_format not in self._areaAverage.keys():
        #         #     self._areaAverage[time_format] = {}
        #         tmp_styles = []
        #         for mask_style,mask in self._masks[grid].items():
        #             tmp_regions = []
        #             for region in mask.region.values:
        #                 tmp_regions.append(np.sum(data * mask.loc[region],axis=(-1,-2)))
        #             tmp_styles.append(xr.concat(tmp_regions, dim='region'))
        #         self._areaAverage[time_format].append(xr.concat(tmp_styles, dim='mask_style'))
        # for time_format in time_formats:
        #     for i in range(2):
        #         self._areaAverage[time_format] = xr.concat(tmp_regions, dim='region')

    def load_area_averages():
        self._tag_combinations = load_pkl(self._working_dir+'meta_info.pkl')
        self._areaAverage = {}
        # for time_format in [ff.split('_')[-2] for ff in glob.glob(self._working_dir+'/areaAverage/*')]:
        for file_name in glob.glob(self._working_dir+'/areaAverage/*'):
            time_format = file_name.split('_')[-2]
            tmp = pd.read_csv(file_name)

            dims = [dim for dim in tmp.columns if dim not in ['Unnamed: 0','value','time']] + ['time']
            coords = {dim:np.unique(tmp[dim])  for dim in dims}
            data = xr.DataArray(data = np.zeros([len(coords[dim]) for dim in dims])*np.nan, coords=coords, dims=dims)

            for grid,tag_combis1 in self._tag_combinations.items():
                for time_format,tag_combis in {key:val for key,val in tag_combis1.items() if key == time_format}.items():
                    for tag_combi in tag_combis.values():
                        slice = tmp.copy()
                        for tag,val in tag_combi.items():
                            slice = slice.loc[(slice[tag]==val)]
                        for mask_style in np.unique(tmp.mask_style):
                            slice_ = slice.loc[(slice['mask_style']==mask_style)]
                            dims = [dim for dim in slice_.columns if dim not in ['Unnamed: 0','value','time']] + ['time']
                            coords = {dim:np.unique(slice_[dim])  for dim in dims}
                            slice_ = slice_['value'].values
                            for i in range(len(dims)-1):
                                slice_ = np.expand_dims(slice_,0)
                            data = xr.DataArray(data = slice_, coords=coords, dims=dims)

                            if time_format not in self._areaAverage.keys():
                                self._areaAverage[time_format] = data
                            else:
                                self._areaAverage[time_format] = xr.align(self._areaAverage[time_format],data, join='outer')[0]
                            print(self._areaAverage[time_format].shape)


    def plot_transient_csv(self,region,var_name,tags):

        region = 'BEN'
        tags = {'scenario':'hist','experiment':'CORDEX','var_name':'mslp'}

        tmp = xr.Dataset({'value':COU._area_averages[region]}).to_dataframe()
        tmp.reset_index(inplace=True)
        for key,val in tags.items():
            tmp = tmp.loc[(tmp[key])==val]
        tmp = tmp.to_xarray()
        plt.plot(tmp.time,tmp.value)
        plt.savefig(COU._working_dir+'plots/test.png')


    def plot_transient(self,region,tags):
        region = 'BEN'
        tags = {'scenario':'hist','experiment':'CORDEX','var_name':'mslp'}

        tmp = COU._area_averages[region]
        indices = []
        for key in tmp.dims[:-1]:
            if key in tags.keys():
                indices.append(tags[key])
            else:
                indices.append(tmp[key].values)
        tmp = tmp.loc[tuple(indices)]

        plt.close()
        # plt.plot(tmp.time.values,tmp.mean(axis=0))
        # plt.fill_between(tmp.time.values,tmp.min(axis=0),tmp.min(axis=0))

        yearly = tmp.groupby('time.year').max('time')
        plt.plot(yearly.year.values,yearly.mean(axis=0))

        plt.savefig(COU._working_dir+'plots/test.png')


        '''
        annual mean
        '''
        tmp.groupby('time.year').max('time')





COU = country_analysis(iso='BEN', working_directory='/Users/peterpfleiderer/Projects/regioClim_2020/cou_data/BEN')
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


COU.harmonize_data()

COU.load_data()
asda
COU.area_average()

asdas

# COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_002_prmsl.nc',var_name='var2',given_var_name='mslp',tag='rcp85')
# COU.zoom_data(input_file='/Users/peterpfleiderer/Projects/data/JRA55/mon_JRA55_vws.nc',var_name='var33',given_var_name='VWS',tag='JRA55')






# COU = country_analysis(iso='BEN', working_directory='/home/pepflei/regioClim_2020/cou_data/BEN')
# # COU.download_shapefile()
# COU.load_shapefile()
# # COU.create_masks(input_file='/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_ECEARTH_rcp45_.nc4', mask_style='pop2015',add_mask_file = '/home/pepflei/CA/masks/population/population_1990_incrLat.nc', add_mask_name='pop2015', add_mask_var='mask')
# # COU.create_masks(input_file='/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_ECEARTH_rcp45_.nc4', mask_style='latWeight', lat_weighted=True)
# COU.load_mask()
#
# # for file_name in glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/pr/mon_pr_*_*_.nc4'):
# #     model = file_name.split('_')[-3]
# #     scenario = file_name.split('_')[-2]
# #     COU.zoom_data(input_file=file_name,var_name='pr',given_var_name='pr',tags={'scenario':scenario,'experiment':'CORDEX','model':model})
# #
# # for file_name in glob.glob('/p/projects/ikiimp/RCM_BC/ISIMIP2b_bc/GCMinput/monthly/tas/mon_tas_*_*_.nc4'):
# #     model = file_name.split('_')[-3]
# #     scenario = file_name.split('_')[-2]
# #     COU.zoom_data(input_file=file_name,var_name='tas',given_var_name='tas',tags={'scenario':scenario,'experiment':'CORDEX','model':model})
#
#
# COU.load_data()
# COU.area_average()
# asdas


#
