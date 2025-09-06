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
    
def get_era5_hourly(dataset_name='reanalysis-era5-land', 
             var=["10m_u_component_of_wind","10m_v_component_of_wind"], 
             year=None,
             month=None,
             day=None,
             time=None,
             grid=[1.0, 1.0],
             area=[90, -180, -90, 180],
             download_flag = False,
             download_file='./output.nc'
            ):
    
    import cdsapi
    import xarray as xr
    from urllib.request import urlopen
    
    # start the cdsapi client
    c = cdsapi.Client()
    if time is None:
        time = [
            "00:00", "01:00", "02:00",
            "03:00", "04:00", "05:00",
            "06:00", "07:00", "08:00",
            "09:00", "10:00", "11:00",
            "12:00", "13:00", "14:00",
            "15:00", "16:00", "17:00",
            "18:00", "19:00", "20:00",
            "21:00", "22:00", "23:00"
        ]

    params = {
    "variable": var,
    "year": year,
    "month": month,
    "day": day,
    "time": time,
    "data_format": "netcdf",
    "download_format": "unarchived",
    "grid": grid,
    "area": area,
    }
        
    # file object
    fl=c.retrieve(dataset_name, params) 
    
    # download the file 
    if download_flag:
        fl.download(f"{download_file}")
    
    # load into memory and return xarray dataset
    with urlopen(fl.location) as f:
        return xr.open_dataset(f.read())
    
from math import radians, sin, cos, sqrt, atan2

# shortest distance between two points on earth
def haversine(lat1, lon1, lat2, lon2):
    R = 6371.0  # Earth radius in km
    phi1, phi2 = radians(lat1), radians(lat2)
    dphi = radians(lat2 - lat1)
    dlambda = radians(lon2 - lon1)
    a = sin(dphi/2)**2 + cos(phi1)*cos(phi2)*sin(dlambda/2)**2
    c = 2 * atan2(sqrt(a), sqrt(1 - a))
    return R * c

import numpy as np
# line between two points on sphere
def interpolate_great_circle(lat1, lon1, lat2, lon2, n_points):
    # Convert to radians
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    # Compute delta
    d = 2 * np.arcsin(np.sqrt(
        np.sin((lat2 - lat1)/2)**2 +
        np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1)/2)**2
    ))
    # Interpolate
    f = np.linspace(0, 1, n_points)
    A = np.sin((1 - f) * d) / np.sin(d)
    B = np.sin(f * d) / np.sin(d)
    x = A * np.cos(lat1) * np.cos(lon1) + B * np.cos(lat2) * np.cos(lon2)
    y = A * np.cos(lat1) * np.sin(lon1) + B * np.cos(lat2) * np.sin(lon2)
    z = A * np.sin(lat1) + B * np.sin(lat2)
    lats = np.arctan2(z, np.sqrt(x**2 + y**2))
    lons = np.arctan2(y, x)
    return np.degrees(lats), np.degrees(lons)

# Calculate bearing (forward azimuth) from start to end point
def calculate_bearing(lat1, lon1, lat2, lon2):
    """
    Calculate the bearing between two points.
    All args in degrees.
    Returns bearing in degrees (0° = North, 90° = East).
    """
    lat1_rad = np.radians(lat1)
    lat2_rad = np.radians(lat2)
    delta_lon = np.radians(lon2 - lon1)
    x = np.sin(delta_lon) * np.cos(lat2_rad) # East-West component
    y = np.cos(lat1_rad) * np.sin(lat2_rad) - np.sin(lat1_rad) * np.cos(lat2_rad) * np.cos(delta_lon) # North-South component
    bearing_rad = np.arctan2(x, y)
    bearing_deg = (np.degrees(bearing_rad) + 360) % 360
    return bearing_deg