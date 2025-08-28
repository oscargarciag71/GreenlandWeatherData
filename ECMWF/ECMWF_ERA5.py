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
    dates=['2020-05-01', '2020-06-01', '2020-07-01', '2021-05-01', '2021-06-01', '2021-07-01', '2022-05-01', '2022-06-01', '2022-07-01', '2023-05-01', '2023-06-01', '2023-07-01'],
    grid=[0.25, 0.25],
    area=[90, -90, 50, -0] # N, W, S, E
)

# Extract data
lon = ds_out['longitude'].values
lat = ds_out['latitude'].values
t2m = ds_out['t2m'].squeeze().values - 273.15  # Kelvin → Celsius
# %%
# --- create figure and axes ---
from matplotlib.colors import BoundaryNorm

# --- Combined figure with two subplots ---
fig, axes = plt.subplots(
    4, 3, figsize=(18, 18),
    subplot_kw={'projection': ccrs.NorthPolarStereo()}
)

# Define boundaries for color levels
boundaries = np.linspace(-25, 15, 41)  # Example: 41 intervals from -25 to 15
norm = BoundaryNorm(boundaries, ncolors=plt.get_cmap('coolwarm').N, clip=True)

for i, ax in enumerate(axes.flat):
    # Add features
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3)
    # ax.set_global() # global map extend
    ax.set_extent([-75, -10, 65, 85], crs=ccrs.PlateCarree())
    # Plot data
    map1 = ax.contourf(
        lon, lat, t2m[i, :, :],
        transform=ccrs.PlateCarree(),
        # levels=30, cmap='coolwarm', vmin=-30, vmax=30, extend='both'
        levels=boundaries, cmap='coolwarm', norm=norm, extend='both'
    )
    # Title
    month = ["May 2020", "June 2020", "July 2020", "May 2021", "June 2021", "July 2021", "May 2022", "June 2022", "July 2022", "May 2023", "June 2023", "July 2023"][i]
    ax.set_title(f'{month}', fontsize=14, y=1.05)

# # Shared colorbar
# cbar = fig.colorbar(map1, ax=axes, orientation='vertical', shrink=1, pad=0.08)

# Shared colorbar using boundaries
cbar = fig.colorbar(map1, ax=axes, orientation='vertical', shrink=1, pad=0.08, boundaries=boundaries, ticks=boundaries[::2])

cbar.set_label('°C', size=12, rotation=0)
cbar.ax.tick_params(labelsize=12)

fig.suptitle("ERA5 Monthly Mean 2m Temperature", fontsize=20, x=0.4, y=0.93)

# plt.tight_layout()
plt.savefig(r"C:\Users\SchwarzN\OneDrive - Université de Fribourg\Private\2026_Greenland\WeatherAnalysis\GreenlandWeatherData\ECMWF\era5_monthly_mean_2m_temperature.png", dpi=300, bbox_inches="tight")
plt.show()

# %%
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Example: central Greenland point
lon_pt, lat_pt = -40.0, 72.0

# Download 10m wind components
ds_wind = get_era5(
    dataset_name="reanalysis-era5-single-levels",
    var=["10m_u_component_of_wind", "10m_v_component_of_wind"],
    dates=["2021-06-01", "2021-06-10"],
    grid=[0.25, 0.25]
)
# %%
# compute wind speed
u10 = ds_wind["u10"].squeeze()
v10 = ds_wind["v10"].squeeze()
wind_speed = np.sqrt(u10**2 + v10**2)

# --------------------
# 1) Time series at central Greenland
# --------------------
# Find nearest gridpoint
pt = ds_wind.sel(latitude=lat_pt, longitude=lon_pt, method="nearest")

u_pt = pt["u10"].values
v_pt = pt["v10"].values
ws_pt = np.sqrt(u_pt**2 + v_pt**2)

plt.figure(figsize=(10,4))
plt.plot(pt["valid_time"], ws_pt, label="Wind speed (m/s)", color="tab:blue")
plt.axhline(5, color="gray", linestyle="--", label="~5 m/s (kiteable)")
plt.legend()
plt.title("10m Wind Speed at Central Greenland (72N, 40W)")
plt.ylabel("m/s")
plt.grid()
plt.show()

# --------------------
# 2) Map plot for June monthly mean
# --------------------
fig, ax = plt.subplots(
    figsize=(7,6),
    subplot_kw={"projection": ccrs.NorthPolarStereo()}
)

ax.set_extent([-75, -10, 59, 85], crs=ccrs.PlateCarree())

# monthly mean wind speed
ws_mean = wind_speed.mean(dim="valid_time")

map1 = ax.contourf(
    ds_wind["longitude"], ds_wind["latitude"], ws_mean,
    transform=ccrs.PlateCarree(),
    cmap="viridis", levels=20
)

plt.colorbar(map1, ax=ax, orientation="horizontal", pad=0.05, label="m/s")
ax.coastlines()
ax.set_title("ERA5 Monthly Mean 10m Wind Speed – June 2021")
plt.show()

# %%
