# -*- coding: utf-8 -*-
'''
Class to analyze climate data on national (& sub national) scale
'''
# Author: Peter Pfleiderer <peter.pfleiderer@climateanalytics.org>
#
# License: GNU General Public License v3.0

import sys,glob,os,pickle,time
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
import threading

os.environ['OPENBLAS_NUM_THREADS'] = '1'

def save_pkl(obj, name ):
	with open(name, 'wb') as f:
		pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_pkl(name ):
	with open( name, 'rb') as f:
		return pickle.load(f)

class shaped_object(object):

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
		if os.path.isdir(self._working_dir+'/raw_data')==False: os.system('mkdir '+self._working_dir+'/raw_data')
		if os.path.isdir(self._working_dir+'/areaAverage')==False:os.system('mkdir '+self._working_dir+'/areaAverage')

		self._masks={}
		self._region_names = {}
		self._region_polygons={}
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

	def read_shapefile(self, shape_file, long_name=None, short_name=None):
		adm_shapefiles=shapereader.Reader(shape_file).records()

		if long_name is None:
			print('Please specify which variables should be used as "long_name" and as "short_name"')
			for item in adm_shapefiles:
				shape,region=item.geometry,item.attributes
				print(item.attributes)
				return 0

		# collect all shapes of region
		for item in adm_shapefiles:
			shape,region=item.geometry,item.attributes
			name_full = u''+region[long_name].replace('_','-')
			name=u''+region[short_name].replace('_','-')
			self._region_names[name]=name_full
			# simplify could be added here to speed up things
			try:
				self._region_polygons[name]=MultiPolygon(shape)
			except:
				self._region_polygons[name]=Polygon(shape)


		save_pkl(self._region_names, self._working_dir+'region_names.pkl')
		save_pkl(self._region_polygons, self._working_dir+'region_polygons.pkl')

	def load_shapefile(self):
		self._region_names = load_pkl(self._working_dir+'region_names.pkl')
		self._region_polygons = load_pkl(self._working_dir+'region_polygons.pkl')


	##########
	# MASKS
	##########

	def load_mask(self):
		for mask_file in glob.glob(self._working_dir+'/masks/*.nc4'):
			small_grid = mask_file.split('_')[-2]
			if small_grid not in self._masks.keys():
				self._masks[small_grid]={}

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
					if mask_style=='latWeight':
						# multiply overlap with latitude weighting
						overlap[j,i] = intersect*np.cos(np.radians(lat[j]))

		# renormalize overlap to get sum(mask)=1
		overlap_sum=np.nansum(overlap.copy().flatten())
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

	def zoom_data(self,input_file,var_name,given_var_name=None,tags={}, lat_name='lat', lon_name='lon', time_cutoff=None, input_array=None):
		lon,lat,grid,lon_shift = self.identify_grid(input_file,lat_name,lon_name)
		if input_array is None:
			_input = xr.open_dataset(input_file)[var_name].squeeze()
		else:

			_input = input_array

		mask = self._masks[self._grid_dict[grid]][list(self._masks[self._grid_dict[grid]].keys())[0]]

		# clean data
		_input.time.values = np.asarray(_input.time,'datetime64[s]')
		med_diff = np.median(np.diff(_input.time,1))
		time_steps = np.array([0] + [i+1 for i,tt in enumerate(np.diff(_input.time,1)) if np.abs(int(tt) - med_diff) < med_diff])
		_input = _input[time_steps,:,:]

		# zoom
		zoomed_data = _input.loc[:,mask.lat,mask.lon]

		# select
		if time_cutoff is not None:
			zoomed_data.sel(time=slice('2006-01-01','2100-01-01'))

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

	def scout_data(self):
		start = time.time()
		# get information about all files
		grids = [nn.split('_')[-2] for nn in glob.glob(self._working_dir+'/masks/*.nc4')]
		self._tags = ['var_name']
		tag_combinations = {grid:{} for grid in grids}
		coords_dict = {grid:{'monthly':{'var_name':[]},'yearly':{'var_name':[]},'daily':{'var_name':[]}} for grid in grids}
		coords_dict_area = {'monthly':{'var_name':[]},'yearly':{'var_name':[]},'daily':{'var_name':[]}}
		for file_name in glob.glob(self._working_dir+'/raw_data/*.nc4'):
			nc_in = xr.open_dataset(file_name)
			var_name = file_name.split('/')[-1].split('_')[2]
			grid = nc_in.attrs['grid']

			# identify and harmonize time format
			tmp_time,time_format,index = self.harmonize_time(nc_in['time'])

			# store coordinates
			if time_format is not None:
				tmp_time = xr.DataArray(tmp_time, coords={'time':tmp_time}, dims=['time'])
				if 'time' not in coords_dict[grid][time_format].keys():
					coords_dict[grid][time_format]['time'] = tmp_time
					coords_dict_area[time_format]['time'] = tmp_time
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

				if time_format not in tag_combinations[grid].keys():
					tag_combinations[grid][time_format] = {}
				tag_combinations[grid][time_format][file_name] = tag_combi

		save_pkl(tag_combinations, self._working_dir+'meta_info.pkl')
		save_pkl(coords_dict, self._working_dir+'grided_coords.pkl')
		for time_format,coords in coords_dict_area.items():
			for key,val in coords.items():
				coords_dict_area[time_format][key] = np.unique(val)
		save_pkl(coords_dict_area, self._working_dir+'coords_areaAverage.pkl')

		print(time.time() - start); start = time.time()

	def read_griddedData(self):
		coords_dict = load_pkl(self._working_dir+'grided_coords.pkl')
		tag_combinations = load_pkl(self._working_dir+'meta_info.pkl')

		# create data arrays
		self._data = {grid:{} for grid in coords_dict.keys()}
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
					self._data[grid][time_format] = xr.DataArray(data=np.zeros([len(coords[key]) for key in dims])*np.nan, coords=coords, dims=dims)

		# fill data arrays
		for grid,tmp1_dict in tag_combinations.items():
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
						self._data[grid][time_format].loc[nc_in.time,nc_in.lat,nc_in.lon] = nc_in[var_name]
					if len(indices) == 1:
						self._data[grid][time_format].loc[indices[0],nc_in.time,nc_in.lat,nc_in.lon] = nc_in[var_name]
					if len(indices) == 2:
						self._data[grid][time_format].loc[indices[0],indices[1],nc_in.time,nc_in.lat,nc_in.lon] = nc_in[var_name]
					if len(indices) == 3:
						self._data[grid][time_format].loc[indices[0],indices[1],indices[2],nc_in.time,nc_in.lat,nc_in.lon] = nc_in[var_name]
					if len(indices) == 4:
						self._data[grid][time_format].loc[indices[0],indices[1],indices[2],indices[3],nc_in.time,nc_in.lat,nc_in.lon] = nc_in[var_name]

				if len(list(tmp2_dict.keys())) > 0:
					self._data[grid][time_format].sortby(['lat', 'lon','time'])

	def print_data(self):
		print('Loaded data into self._data')
		for grid, tmp1_dict in self._data.items():
			for time_format, tmp in tmp1_dict.items():
				print('grid: '+grid)
				print('time format: '+time_format)
				print(tmp.coords,'\n')
				xr.Dataset({'data':tmp}).to_netcdf(self._working_dir+'gridded_data/'+'_'.join([self._iso,grid,time_format])+'.nc')

	def calculate_areaAverage(self, regions=None):
		tag_combinations = load_pkl(self._working_dir+'meta_info.pkl')
		coords_area = load_pkl(self._working_dir+'coords_areaAverage.pkl')
		time_formats = np.array([[time_format for time_format in tmp1.keys()] for tmp1 in tag_combinations.values()]).flatten()

		if regions is None:
			regions = self._region_names.keys()

		for region in regions:
			for time_format in time_formats:
				dims = [dim for dim in coords_area[time_format].keys() if dim not in ['lat','lon','time']] + ['time']
				coords = {dim:coords_area[time_format][dim]  for dim in dims}
				coords['mask_style'] = np.array([list(tmp1.keys()) for grid,tmp1 in self._masks.items()]).flatten()
				dims = ['mask_style'] + dims
				tmp_xr = xr.DataArray(data=np.zeros([len(coords[key]) for key in dims])*np.nan, coords=coords, dims=dims)

				for grid,tmp1 in tag_combinations.items():
					for time_format,tmp2 in tmp1.items():
						for tag_combi in tmp2.values():
							for mask_style,mask in self._masks[grid].items():
								# if tag_combi['var_name'] == 'tas':
								# 	asdas

								tmp = self._data[grid][time_format].loc[tuple(tag_combi[dim] for dim in self._data[grid][time_format].dims[:-3])].copy()
								nans = np.isnan(tmp)

								tmp *= mask.loc[region]
								av = tmp.sum(axis=(-2,-1))
								av.values[np.all(nans, axis=(1,2))] = np.nan

								tag_combi_ = tag_combi.copy()
								tag_combi_['region'] = region
								tag_combi_['mask_style'] = mask_style
								indices = [tag_combi_[dim] for dim  in tmp_xr.dims if dim not in ['time']]

								if len(indices) == 0:
									tmp_xr.loc[av.time] = av
								if len(indices) == 1:
									tmp_xr.loc[indices[0],av.time] = av
								if len(indices) == 2:
									tmp_xr.loc[indices[0],indices[1],av.time] = av
								if len(indices) == 3:
									tmp_xr.loc[indices[0],indices[1],indices[2],av.time] = av
								if len(indices) == 4:
									tmp_xr.loc[indices[0],indices[1],indices[2],indices[3],av.time] = av
								if len(indices) == 5:
									tmp_xr.loc[indices[0],indices[1],indices[2],indices[3],indices[4],av.time] = av
								if len(indices) == 6:
									tmp_xr.loc[indices[0],indices[1],indices[2],indices[3],indices[4],indices[5],av.time] = av

				if os.path.isfile(self._working_dir+'areaAverage/'+region+'_'+time_format+'_areaAverage.nc'):
					os.system('rm '+self._working_dir+'areaAverage/'+region+'_'+time_format+'_areaAverage.nc')
				xr.Dataset({'areaAverage':tmp_xr}).to_netcdf(self._working_dir+'areaAverage/'+region+'_'+time_format+'_areaAverage.nc')

	def load_areaAverage(self):
		tag_combinations = load_pkl(self._working_dir+'meta_info.pkl')
		self._areaAverage = {}
		time_formats = np.array([[time_format for time_format in tmp1.keys()] for tmp1 in tag_combinations.values()]).flatten()
		self._areaAverage = {}
		for time_format in time_formats:
			tmp = xr.open_mfdataset(self._working_dir+'areaAverage/*_'+time_format+'_areaAverage.nc', concat_dim='region', combine='nested')
			tmp.coords['region'] = [ff.split('/')[-1].split('_')[0] for ff in glob.glob(self._working_dir+'areaAverage/*_'+time_format+'_areaAverage.nc')]
			self._areaAverage[time_format] = tmp['areaAverage'].load()

	def select_areaAverage(self, time_format, tags):
		tmp = self._areaAverage[time_format]
		indices = []
		for key in tmp.dims[:-1]:
			if key in tags.keys():
				indices.append(tags[key])
			else:
				indices.append(tmp[key].values)
		return tmp.loc[tuple(indices)]

	def select_griddedData(self, time_format, tags):
		indices = []
		for grid,tmp in self._data.items():
			found_dims = []
			for key in tags.keys():
				if key in tmp[time_format].dims:
					found_dims.append(key)
			if len(found_dims) == len(list(tags.keys())):
				for key in tmp[time_format].dims[:-2]:
					if key in tags.keys():
						indices.append(tags[key])
					else:
						indices.append(tmp[time_format][key].values)
				return tmp[time_format].loc[tuple(indices)]

	def unit_conversion(self, time_format, var_name, addition=0, factor=1):
		tags = {'var_name':var_name}
		indices = []
		for grid,tmp in self._data.items():
			found_dims = []
			for key in tags.keys():
				if key in tmp[time_format].dims:
					found_dims.append(key)
			if len(found_dims) == len(list(tags.keys())):
				for key in tmp[time_format].dims[:-2]:
					if key in tags.keys():
						indices.append(tags[key])
					else:
						indices.append(tmp[time_format][key].values)

				print(indices)
				print(np.nanmedian(tmp[time_format].loc[tuple(indices)]))
				tmp[time_format].loc[tuple(indices)] += addition
				tmp[time_format].loc[tuple(indices)] *= factor
				print(np.nanmedian(tmp[time_format].loc[tuple(indices)]))






#
