# %%
import cdsapi
from urllib.request import urlopen

c = cdsapi.Client()

wind_box = c.retrieve(
    'reanalysis-era5-land',
    {
        'variable': ['10m_u_component_of_wind', '10m_v_component_of_wind'],
        'year': ['2020', '2021'],
        'month': ['01', '02', '03'],
        'day': ['01', '02', '03'],
        'time': ['00:00', '06:00', '12:00', '18:00'],
        'area': [85, -75, 59, -10],  # N, W, S, E (Greenland box)
        'format': 'netcdf',
    },
    'greenland_wind.nc')
wind_box.download(f"{'./output.nc'}")
with urlopen(wind_box.location) as f:
    ds = xr.open_dataset(f.read())
# %%
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Load the NetCDF file
ds = xr.open_dataset("greenland_wind.nc")

# Variables
u = ds['u10']  # eastward wind
v = ds['v10']  # northward wind

# Select one timestep (first one here)
u0 = u.isel(time=0)
v0 = v.isel(time=0)

# Plot wind vectors on a map
plt.figure(figsize=(8, 10))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.coastlines()

# Subsample grid for clarity
step = 4
plt.quiver(
    u0['longitude'][::step], u0['latitude'][::step],
    u0.values[::step, ::step], v0.values[::step, ::step],
    transform=ccrs.PlateCarree(), scale=400
)

plt.title(f"Wind at {str(u0['time'].values)} (10m)")
plt.show()
