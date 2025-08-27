# %%
def get_era5(dataset_name='reanalysis-era5-single-levels', 
             var=None, 
             dates=None,
             pressure_level=None,
             grid=[1.0, 1.0],
             area=[90, -180, -90, 180],
             download_flag = False,
             download_file='./output.nc'
            ):
    ''' Get ERA5 reanalysis output from the web
    this script grabs ERA5 variables from the web and stores them 
    in an xarray dataset. 
    
    the ERA5 CDS API must be installed on the local machine.
    See section 4 here: https://cds.climate.copernicus.eu/api-how-to
    
    Parameters
    ----------                  
    dataset_name: str, default 'reanalysis-era5-single-levels'
        name of dataset to use. Options include:
        * 'reanalysis-era5-single-levels'
        * 'reanalysis-era5-single-levels-monthly-means'
        * 'reanalysis-era5-pressure-levels'
        * 'reanalysis-era5-pressure-levels-monthly-means'
        * 'reanalysis-era5-land'
        * 'reanalysis-era5-land-monthly-means'
        
    dates: list of strings or datetime64, default None
        example ['1980-01-01', '2020-12-31']
    var: str, default None
        name of variable to download
        example '2m_temperature'
    pressure_level: str, default None
        pressure level to grab data on
    grid: list, deafult [1.0, 1.0]
        spatial lat, lon grid resolution in deg
    area: list, default [90,-180,-90, 180]
        area extent download [N, W, S, E]
    download_flag = True or False, default False
        flag to download data or not
    download_file= str, default './output.nc'
        path to where data should be downloaed to.
        data only downloaded if download_flag is True
    Returns
    -------
    ds: xarrayDataSet
        all the data will be in an xarray dataset
        
    Example
    -------
    ds = get_era5(dataset_name='reanalysis-era5-single-levels-monthly-means', 
                 var='2m_temperature', 
                 dates=['2021-02-01'],
                 grid=[0.25, 0.25])
        
    Notes
    -------    
    # cdsapi code is here
    https://github.com/ecmwf/cdsapi/tree/master/cdsapi
    # information on api is here
    https://confluence.ecmwf.int/display/CKB/Climate+Data+Store+%28CDS%29+API+Keywords
    # era5 dataset information is here
    https://confluence.ecmwf.int/display/CKB/The+family+of+ERA5+datasets
    '''
    import cdsapi
    import xarray as xr
    import pandas as pd
    from urllib.request import urlopen
    
    # test if acceptable pressure level
    acceptable_pressures = [1, 2, 3, 5, 7, 10, 20, 30, 50, 70, range(100, 1000, 25)]
    if pressure_level not in [str(lev) for lev in acceptable_pressures]:
        print(f"!! Pressure level must be in this list: {acceptable_pressures}")
    
    # start the cdsapi client
    c = cdsapi.Client()

    # parameters
    params = dict(
        format = "netcdf",
        product_type = "reanalysis",
        variable = var,
        grid = grid,
        area = area,
        date = list(dates.strftime('%Y-%m-%d %H:%M')) \
               if isinstance(dates, pd.core.indexes.datetimes.DatetimeIndex)\
               else dates,
        )

    # what to do if asking for monthly means
    if dataset_name in ["reanalysis-era5-single-levels-monthly-means", 
                        "reanalysis-era5-pressure-levels-monthly-means",
                        "reanalysis-era5-land-monthly-means"]:
        params["product_type"] = "monthly_averaged_reanalysis"
        _ = params.pop("date")
        params["time"] = "00:00"
        
        # if time is in list of pandas format
        if isinstance(dates, list):
            dates_pd = pd.to_datetime(dates)
            params["year"] = sorted(list(set(dates_pd.strftime("%Y"))))
            params["month"] = sorted(list(set(dates_pd.strftime("%m"))))
        else:
            params["year"] = sorted(list(set(dates.strftime("%Y"))))
            params["month"] = sorted(list(set(dates.strftime("%m"))))
            
        
    # if pressure surface
    if dataset_name in ["reanalysis-era5-pressure-levels-monthly-means",
                        "reanalysis-era5-pressure-levels"]:
        params["pressure_level"] = pressure_level
    
    # product_type not needed for era5_land
    if dataset_name in ["reanalysis-era5-land"]:
        _ = params.pop("product_type")
        
    # file object
    fl=c.retrieve(dataset_name, params) 
    
    # download the file 
    if download_flag:
        fl.download(f"{download_file}")
    
    # load into memory and return xarray dataset
    with urlopen(fl.location) as f:
        return xr.open_dataset(f.read())
# %%
'''
this is an example os how to use the get_era5() function
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcdefaults()  
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# --- get the data ---
# Requires you to have get_era5() available
# https://gist.github.com/lgloege/f461f8d192e99fe7c36760a7a856b007
ds_out = get_era5(
    dataset_name='reanalysis-era5-single-levels-monthly-means',
    var='2m_temperature',
    dates=['2021-02-01'],
    grid=[0.25, 0.25]
)

# Extract data
lon = ds_out['longitude'].values
lat = ds_out['latitude'].values
t2m = ds_out['t2m'].squeeze().values - 273.15  # Kelvin → Celsius

# --- create figure and axes ---
fig, ax = plt.subplots(
    figsize=(10, 5),
    subplot_kw={'projection': ccrs.Robinson(central_longitude=0)}
)

# add features
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linewidth=0.3)
ax.set_global()

# plot data
map1 = ax.contourf(
    lon, lat, t2m,
    transform=ccrs.PlateCarree(),
    levels=20, cmap='coolwarm', vmin=-30, vmax=30, extend='both'
)

# title
ax.set_title('ERA5 monthly mean 2m temperature - February 2021', fontsize=14)

# colorbar
cbar = plt.colorbar(map1, ax=ax, orientation='horizontal', shrink=0.6, pad=0.05)
cbar.set_label('°C', size=12)
cbar.ax.tick_params(labelsize=12)

plt.show()
