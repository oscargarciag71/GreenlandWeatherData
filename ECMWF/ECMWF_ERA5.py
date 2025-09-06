# %%
import sys
sys.path.append('C:/Users/SchwarzN/OneDrive - Université de Fribourg/Private/2026_Greenland/WeatherAnalysis/GreenlandWeatherData/ECMWF')
import importlib
from utils import get_era5, get_era5_hourly
importlib.reload(sys.modules['utils'])
# %%
'''
Monthly mean temperatures plotted over Greenland
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

# %% Tasiilaq until Ilulissat
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import sys
sys.path.append('C:/Users/SchwarzN/OneDrive - Université de Fribourg/Private/2026_Greenland/WeatherAnalysis/GreenlandWeatherData/ECMWF')
import importlib
from utils import haversine, interpolate_great_circle, calculate_bearing
importlib.reload(sys.modules['utils'])
import matplotlib

# Define start and end points (lat, lon)
start_lat, start_lon = 66.071024, -38.332835 # Tasiilaq
end_lat, end_lon = 69.213215, -49.422260 # Ilulissat

distance_km = haversine(start_lat, start_lon, end_lat, end_lon)

bearing_deg = calculate_bearing(start_lat, start_lon, end_lat, end_lon)
print(f"Bearing from start to end point: {bearing_deg:.2f}°")

# Starting date
start_date = "05-30" # May 30
end_date = "07-08" # July 8
year = 2020
days = np.arange(np.datetime64(f"{year}-{start_date}"), np.datetime64(f"{year}-{end_date}") + np.timedelta64(1, 'D'))
dates = [str(day) for day in days]
n_days = len(dates)
# Number of interpolation points (including endpoints)
n_points = n_days
lats, lons = interpolate_great_circle(start_lat, start_lon, end_lat, end_lon, n_points)

ds_wind = get_era5(
    dataset_name="reanalysis-era5-single-levels", # "reanalysis-era5-land"
    var=["10m_u_component_of_wind", "10m_v_component_of_wind"],
    dates=dates,
    grid=[0.1, 0.1],
    area=[70, -50, 66, -38]
)

# compute wind speed
u10 = ds_wind["u10"].squeeze()
v10 = ds_wind["v10"].squeeze()
wind_speed = np.sqrt(u10**2 + v10**2)
# %%
# --------------------
#
# --------------------
# Find nearest gridpoints

# For each day, select the wind speed at the nearest grid point to the interpolated path
ws_along_path = []
wind_dir_along_path = []
for i in range(n_points):
    u_pt = u10.sel(valid_time=dates[i], latitude=lats[i], longitude=lons[i], method="nearest").values
    v_pt = v10.sel(valid_time=dates[i], latitude=lats[i], longitude=lons[i], method="nearest").values
    ws_pt = np.sqrt(u_pt**2 + v_pt**2)
    # Wind direction in degrees: 0 = North, 90 = East, 180 = South, 270 = West
    wd_pt = (np.arctan2(u_pt, v_pt) * 180 / np.pi) % 360
    ws_along_path.append(ws_pt)
    wind_dir_along_path.append(wd_pt)
ws_along_path = np.array(ws_along_path)
wind_dir_along_path = np.array(wind_dir_along_path)
# Define 24 colors for 24 wind direction bins (15° each)

n_bins = 24
cmap = plt.get_cmap('hsv', n_bins)
angle_bins = np.arange(0, 360, 15)
color_array = []

for angle in wind_dir_along_path:
    bin_idx = int(angle // 15) % n_bins
    rgba = cmap(bin_idx)
    hex_color = matplotlib.colors.rgb2hex(rgba)
    color_array.append(hex_color)

# Assign color to bearing_deg angle
bearing_bin_idx = int(bearing_deg // 15) % n_bins
bearing_color = matplotlib.colors.rgb2hex(cmap(bearing_bin_idx))

color_array = np.array(color_array)

# Plot the wind along the path
plt.figure(figsize=(12, 5))
plt.plot(dates, ws_along_path, label="Wind speed (m/s)", color="tab:blue")
plt.axhline(5, color="gray", linestyle="--", label="~5 m/s (kiteable)")

# Overlay wind direction as arrows
skip = max(1, len(dates)//25)  # Show at most ~25 arrows for clarity
for i in range(0, len(dates), skip):
    # Arrow parameters
    y = ws_along_path[i]
    angle = wind_dir_along_path[i]
    arrow_length = 0.8
    dx = arrow_length * np.sin(np.deg2rad(angle))
    dy = arrow_length * np.cos(np.deg2rad(angle))
    plt.arrow(i, y, dx, dy, head_width=0.5, head_length=0.3, fc=color_array[i], ec=color_array[i],
              length_includes_head=True, alpha=0.7, zorder=3,
              overhang=0.3, linewidth=2)
plt.arrow(13, 8, 3*np.sin(np.deg2rad(bearing_deg)), 3*np.cos(np.deg2rad(bearing_deg)),
          head_width=0.5, head_length=0.3, fc=bearing_color, ec=bearing_color,
          length_includes_head=True, alpha=0.7, zorder=3,
          overhang=0.3, linewidth=4)

plt.legend()
plt.title("10m Wind Speed and Direction Along the Path (Tasiilaq to Ilulissat)")
plt.ylabel("Wind speed (m/s)")
plt.xticks(ticks=range(0, len(dates), skip), labels=[dates[i][5:] for i in range(0, len(dates), skip)], rotation=45)
plt.grid()
plt.tight_layout()
plt.show()

# Plot the interpolation path
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize=(8, 6))
ax.plot(lons, lats, marker='o', color='blue', markersize=2, label='Interpolated Path')
ax.scatter([start_lon], [start_lat], color='green', s=50, label='Start Point')
ax.scatter([end_lon], [end_lat], color='red', s=50, label='End Point')
ax.coastlines()
ax.set_extent([-55, -35, 65, 70])
ax.legend()
ax.set_title('Interpolated Path from Start to End Point')
plt.show()

# %%
#
# Hourly wind analysis along route
#
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import sys
sys.path.append('C:/Users/SchwarzN/OneDrive - Université de Fribourg/Private/2026_Greenland/WeatherAnalysis/GreenlandWeatherData/ECMWF')
import importlib
from utils import get_era5, get_era5_hourly, haversine, interpolate_great_circle, calculate_bearing
importlib.reload(sys.modules['utils'])
import matplotlib

# Define start and end points (lat, lon)
start_lat, start_lon = 69.753550, -50.279133 # Ilulissat
end_lat, end_lon = 65.9654610324651, -38.322433959692717 # Tasiilaq

distance_km = haversine(start_lat, start_lon, end_lat, end_lon)

bearing_deg = calculate_bearing(start_lat, start_lon, end_lat, end_lon)
print(f"Bearing from start to end point: {bearing_deg:.2f}°")

# Starting date
start_date = "05-30" # May 30
end_date = "07-08" # July 8

year = 0000
days = np.arange(np.datetime64(f"{year}-{start_date}"), np.datetime64(f"{year}-{end_date}") + np.timedelta64(1, 'D'))
dates = [str(day) for day in days]
n_days = len(dates)

years = np.arange(1986, 2026)

# Number of interpolation points (including endpoints) set to number of days
n_points = n_days
lats, lons = interpolate_great_circle(start_lat, start_lon, end_lat, end_lon, n_points)

# %%
# Collect wind for all years, dates and hours
hourly_wind_results = {}

for year in years:
    hourly_wind_results[year] = {}
    for date_str in dates:
        # Extract month and day from date string (format: 'YYYY-MM-DD')
        _, month, day = date_str.split('-')
        # Download hourly wind data for this date
        ds_hourly_wind = get_era5_hourly(
            dataset_name="reanalysis-era5-single-levels",
            var=["10m_u_component_of_wind", "10m_v_component_of_wind"],
            year=str(year),
            month=month,
            day=day,
            grid=[0.1, 0.1],
            area=[70, -50, 66, -38]
        )
        # Store result for later analysis
        # Replace the year in date_str with the current year for consistency
        date_parts = date_str.split('-')
        date_with_correct_year = f"{year}-{date_parts[1]}-{date_parts[2]}"
        hourly_wind_results[year][date_with_correct_year] = ds_hourly_wind

# Save each dataset separated by year and date for easy processing afterwards !DON'T SAVE IN ONE-DRIVE DIRECTORY!
for year in years:
    for i in np.arange(len(dates)):
        date = f"{year}-{dates[i].split('-')[1]}-{dates[i].split('-')[2]}"
        ds = hourly_wind_results[year][date]
        ds.to_netcdf(f"C:/Users/SchwarzN/wind_data/hourly_wind_{year}_{date}.nc")
# %% Process data to extract wind speed and direction along the path
import xarray as xr
import numpy as np

# 2D arrays for wind speed and direction: shape (len(years), n_points, 24)
hours = [f"{h:02d}" for h in range(24)]
wind_speed_along_path = np.zeros((len(years), n_points, len(hours)))
wind_dir_along_path = np.zeros((len(years), n_points, len(hours)))

for y_idx, year in enumerate(years):
    for i in np.arange(len(dates)):
        date = f"{year}-{dates[i].split('-')[1]}-{dates[i].split('-')[2]}"
        path = f"C:/Users/SchwarzN/wind_data/hourly_wind_{year}_{date}.nc"
        ds = xr.open_dataset(path)
        u10 = ds["u10"].squeeze()
        v10 = ds["v10"].squeeze()
        # Retrieve for all hours
        for h in range(24):
            u_pt = u10.sel(valid_time=date+"T"+hours[h], latitude=lats[i], longitude=lons[i], method="nearest").values
            v_pt = v10.sel(valid_time=date+"T"+hours[h], latitude=lats[i], longitude=lons[i], method="nearest").values
            ws_pt = np.sqrt(u_pt**2 + v_pt**2)
            # Wind direction in degrees: 0 = North, 90 = East, 180 = South, 270 = West
            wd_pt = (np.arctan2(u_pt, v_pt) * 180 / np.pi) % 360
            wind_speed_along_path[y_idx, i, h] = ws_pt # shape (len(years), len(dates), 24)
            wind_dir_along_path[y_idx, i, h] = wd_pt # shape (len(years), len(dates), 24)

# %% Generate windrose for all locations/date taking into account all years and hours
# ? Is it reasonable to analyze for differences and trends between different years?
import windrose

def plot_windrose(wind_speed, wind_dir, years, date, lat, lon):
    date = f"{date.split('-')[2]}.{date.split('-')[1]}."
    ax = windrose.WindroseAxes.from_ax()
    ax.bar(wind_dir.flatten(), wind_speed.flatten(), normed=True, opening=0.8, edgecolor='white')
    ax.set_title(f"Wind Rose for {years[0]} - {years[-1]} ({date})\nLocation: {lat:.2f}°, {lon:.2f}°")
    ax.set_xlabel("Wind Direction (°)")
    ax.set_ylabel("Wind Speed (m/s)")
    
    ax.legend(title="Wind Speed (m/s)", loc='lower left', bbox_to_anchor=(1.1, 0))
    plt.show()

for day_idx, day in enumerate(dates):
    plot_windrose(wind_speed_along_path[:, day_idx, :], wind_dir_along_path[:, day_idx, :], years, dates[day_idx], lats[day_idx], lons[day_idx])

# %%
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from windrose import WindroseAxes, plot_windrose

wind_data = pd.DataFrame(
    {
        "ws": np.random.random(1200) * 6,
        "wd": np.random.random(1200) * 360,
        "month": np.repeat(range(1, 13), 100),
    }
)


def plot_windrose_subplots(data, *, direction, var, color=None, **kwargs):
    """wrapper function to create subplots per axis"""
    ax = plt.gca()
    ax = WindroseAxes.from_ax(ax=ax)
    plot_windrose(direction_or_df=data[direction], var=data[var], ax=ax, **kwargs)


# this creates the raw subplot structure with a subplot per value in month.
g = sns.FacetGrid(
    data=wind_data,
    # the column name for each level a subplot should be created
    col="month",
    # place a maximum of 3 plots per row
    col_wrap=3,
    subplot_kws={"projection": "windrose"},
    sharex=False,
    sharey=False,
    despine=False,
    height=3.5,
)

g.map_dataframe(
    plot_windrose_subplots,
    direction="wd",
    var="ws",
    normed=True,
    # manually set bins, so they match for each subplot
    bins=(0.1, 1, 2, 3, 4, 5),
    calm_limit=0.1,
    kind="bar",
)

# make the subplots easier to compare, by having the same y-axis range
y_ticks = range(0, 17, 4)
for ax in g.axes:
    ax.set_legend(
        title=r"$m \cdot s^{-1}$", bbox_to_anchor=(1.15, -0.1), loc="lower right"
    )
    ax.set_rgrids(y_ticks, y_ticks)

# adjust the spacing between the subplots to have sufficient space between plots
plt.subplots_adjust(wspace=-0.2)
# %%
# --- Windrose subplots for each date using calculated data ---
import matplotlib.pyplot as plt
from windrose import WindroseAxes

num_dates = len(dates)
cols = 3
rows = int(np.ceil(num_dates / cols))
fig, axes = plt.subplots(rows, cols, subplot_kw={"projection": "windrose"}, figsize=(cols*4, rows*4))
axes = axes.flatten()

for day_idx, day in enumerate(dates):
    ax = axes[day_idx]
    ws = wind_speed_along_path[:, day_idx, :].flatten()
    wd = wind_dir_along_path[:, day_idx, :].flatten()
    ax.bar(wd, ws, normed=True, opening=0.8, edgecolor='white')
    ax.set_title(f"{day} ({lats[day_idx]:.2f}°, {lons[day_idx]:.2f}°)")
    ax.set_xlabel("Wind Direction (°)")
    ax.set_ylabel("Wind Speed (m/s)")
    ax.legend(title="Wind Speed (m/s)", loc='lower left', bbox_to_anchor=(1.1, 0))

# Hide unused subplots if any
for i in range(num_dates, len(axes)):
    fig.delaxes(axes[i])

plt.tight_layout()
plt.show()
# %%
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from windrose import WindroseAxes, plot_windrose

bearing_angle = bearing_deg  # arrow in windrose plots

# --- Prepare your data into a flat dataframe for seaborn ---
records = []
for day_idx, date in enumerate(dates):
    ws = wind_speed_along_path[:, day_idx, :].flatten()
    wd = wind_dir_along_path[:, day_idx, :].flatten()
    lat = lats[day_idx]
    lon = lons[day_idx]
    
    for s, d in zip(ws, wd):
        records.append({
            "ws": s,
            "wd": d,
            "date": date,
            "lat": lat,
            "lon": lon,
        })

wind_data = pd.DataFrame.from_records(records)

# --- Windrose plotting function ---
def plot_windrose_subplots(data, *, direction, var, **kwargs):
    ax = plt.gca()
    ax = WindroseAxes.from_ax(ax=ax)
    plot_windrose(
        direction_or_df=data[direction],
        var=data[var],
        ax=ax,
        **kwargs
    )
    return ax

# --- Build FacetGrid: one subplot per date ---
g = sns.FacetGrid(
    data=wind_data,
    col="date",
    col_wrap=4,        # adjust if too many columns
    subplot_kws={"projection": "windrose"},
    sharex=False,
    sharey=False,
    despine=False,
    height=3.5,
)

g.map_dataframe(
    plot_windrose_subplots,
    direction="wd",
    var="ws",
    normed=True,
    bins=(0.1, 1, 2, 3, 4, 5),
    calm_limit=0.1,
    kind="bar",
)

# --- Format each subplot ---
y_ticks = range(0, 17, 4)
for ax, (date, subdata) in zip(g.axes.flat, wind_data.groupby("date")):
    lat = subdata["lat"].iloc[0]
    lon = subdata["lon"].iloc[0]
    ax.set_title(f"{date}\n({lat:.2f}°, {lon:.2f}°)")

    # Legend
    ax.set_legend(
        title=r"$m \cdot s^{-1}$", bbox_to_anchor=(1.15, -0.1), loc="lower right"
    )
    ax.set_rgrids(y_ticks, y_ticks)

    # --- Add arrow for bearing angle ---
    theta = np.deg2rad(bearing_angle+180)
    rmax = ax.get_ylim()[1]  # maximum radius
    ax.annotate(
        "", 
        xy=(theta, rmax*0.9),  # arrow tip
        xytext=(theta, 0),     # arrow start
        arrowprops=dict(facecolor="red", edgecolor="red", width=2, headwidth=8)
    )

plt.subplots_adjust(wspace=-0.2, hspace=0.4)
plt.show()
