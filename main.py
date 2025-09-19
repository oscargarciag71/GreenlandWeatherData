import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import utils
import cdsapi
import xarray as xr
import pandas as pd
from urllib.request import urlopen

# ----Locations----
# Station 446500
latitude_446500 = 66.48
longitude_446500 = -46.30
# Station 447500
latitude_447500 = 65.18
longitude_447500 = -43.83
# Waypoint 1
WP1_lat = 69.305186
WP1_lon = -48.060573
# Waypoint 2
WP2_lat = 68.668851
WP2_lon = -45.760206
# Waypoint 3
WP3_lat = 67.905773
WP3_lon = -43.650136
# Waypoint 4
WP4_lat = 67.053953
WP4_lon = -41.173379
# Waypoint 5
WP5_lat = 66.369434
WP5_lon = -39.460040
# -----------------

# Uncomment only the parts of code that are necessary

# Import CSV file
# df = pd.read_csv("data/446500.csv", sep=";")

# Plot histograms
# utils.plot_histograms(df)

# Plot wind rose
# utils.plot_wind_rose(df, 304, 354, 254)


# Download data from ERA5 using Copernicus API
# utils.download_era5_data_wind(WP3_lat, WP3_lon, 2008, 2018, "data/ERA5_WP3")

# Convert the ERA5 data into CSV
# utils.convert_era5_data_wind("data/ERA5_WP3", "data/ERA5_WP3/wind_WP3.csv")

# Import CSV file
# df = pd.read_csv("data/ERA5_WP1/wind_WP1.csv", sep=";")

# Plot wind rose
# utils.plot_wind_rose(df, 304)

# Download data from ERA5 using Copernicus API
# utils.download_era5_data_snow(WP1_lat, WP1_lon, 2022, 2024, "data/SNOW_TEST")

# utils.convert_era5_data_snow("data/SNOW_TEST", "data/SNOW_TEST/test .csv")

df = pd.read_csv("data/SNOW_TEST/test .csv", sep=";")

utils.linear_plot(df)
