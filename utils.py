import pandas as pd
import matplotlib.pyplot as plt
from windrose import WindroseAxes
import numpy as np
import cdsapi
import xarray as xr
import os
import glob


def linear_plot(df):
    plt.figure(figsize=(10, 5))
    plt.plot(df.index, df["skt"], label="skt (ºC)", color="blue")
    plt.xlabel("Time")
    plt.ylabel("Values")
    plt.title("Skt")
    plt.legend()
    plt.grid(True)
    plt.show()


def plot_histograms(df):
    # Filter only May, June, July
    months = {5: "May", 6: "June", 7: "July"}

    for month_num, month_name in months.items():
        subset = df[df["Month"] == month_num]
        print(subset)

        # Create figure with 2 histograms (101 and 301)
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle(f"Histograms for {month_name}", fontsize=16)

        # Histogram for column 101
        temp_data = subset["101"].dropna()
        mean_temp = temp_data.mean()
        std_temp = temp_data.std()
        min_temp = temp_data.min()  # Absolute minimum temperature
        temp_data.hist(bins=30, ax=axes[0], edgecolor="black", range=(-35, 20))
        axes[0].set_title("Mean air temperature")
        axes[0].set_xlabel("ºC")
        axes[0].text(
            0.95,
            0.95,
            f"Mean: {mean_temp:.2f} ºC \nStd: {std_temp:.2f} ºC \n Abs Min: {min_temp:.2f} ºC",
            transform=axes[0].transAxes,
            fontsize=12,
            verticalalignment="top",
            horizontalalignment="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7),
        )

        # axes[0].set_ylabel("Frequency")

        # Histogram for column 301
        wind_data = subset["301"].dropna()
        mean_wind = wind_data.mean()
        std_wind = wind_data.std()
        max_wind = wind_data.max()  # Absolute maximum wind speed
        wind_data.hist(bins=30, ax=axes[1], edgecolor="black", range=(0, 25))
        axes[1].set_title("Mean wind speed")
        axes[1].set_xlabel("m/s")
        mean_wind = wind_data.mean()
        std_wind = wind_data.std()
        # axes[1].set_ylabel("Frequency")
        axes[1].text(
            0.95,
            0.95,
            f"Mean: {mean_wind:.2f} m/s \nStd: {std_wind:.2f} m/s \n Abs Max: {max_wind:.2f} m/s",
            transform=axes[1].transAxes,
            fontsize=12,
            verticalalignment="top",
            horizontalalignment="right",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.7),
        )

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.show()


def plot_wind_rose(df, heading):
    kite_angle_1 = (heading - 55) % 360
    kite_angle_2 = (heading + 55) % 360

    # Months of interest
    months = {5: "May", 6: "June", 7: "July"}

    fig = plt.figure(figsize=(18, 6))

    # Create 3 windrose axes manually
    for i, (month, name) in enumerate(months.items()):
        df_month = df[df["Month"] == month][["301", "365"]].dropna()
        ws = df_month["301"]
        wd = df_month["365"]

        # Define axes rectangle: left, bottom, width, height
        rect = [0.05 + i * 0.32, 0.1, 0.3, 0.8]  # adjust spacing between subplots
        ax = WindroseAxes(fig, rect)
        fig.add_axes(ax)

        ax.bar(
            wd,
            ws,
            opening=0.8,
            edgecolor="white",
            bins=[0, 2, 4, 6, 8, 10],
            normed=True,
        )
        ax.set_title(f"Average wind {name}", fontsize=14)

        # Add radial marker at 124°

        theta = np.deg2rad(-heading + 90)
        ax.plot([theta, theta], [0, ax.get_rmax()], "r--", lw=2)

        theta_1 = np.deg2rad(-kite_angle_1 + 90)
        ax.plot([theta_1, theta_1], [0, ax.get_rmax()], "b--", lw=2)

        theta_2 = np.deg2rad(-kite_angle_2 + 90)
        ax.plot([theta_2, theta_2], [0, ax.get_rmax()], "b--", lw=2)

        # ---- Compute percentage ----
        mask = (ws > 2.5) & ((wd < kite_angle_1) | (wd > kite_angle_2))
        percentage = 100 * mask.sum() / len(ws)

        # Add text annotation inside plot
        ax.text(
            0.5,
            -0.15,
            f"Kiting conditions {percentage:.1f}%",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=12,
            color="black",
        )
    # Add a single legend for all subplots
    ax.set_legend(title="Wind speed (m/s)")

    plt.show()


def download_era5_data_wind(
    latitude, longitude, year_start, year_end, output_folder="data/TEST_ERA5"
):
    c = cdsapi.Client()

    # create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    years = [str(y) for y in range(year_start, year_end + 1)]
    months = ["05", "06", "07"]  # May, June, July

    for year in years:
        for month in months:
            filename = f"wind_{year}_{month}.nc"
            filepath = os.path.join(
                output_folder, filename
            )  # store inside chosen folde
            print(f"Downloading {filename}...")

            c.retrieve(
                "reanalysis-era5-single-levels",
                {
                    "product_type": "reanalysis",
                    "variable": ["10m_u_component_of_wind", "10m_v_component_of_wind"],
                    "year": year,
                    "month": month,
                    "day": [str(d).zfill(2) for d in range(1, 32)],  # all days
                    "time": ["00:00", "06:00", "12:00", "18:00"],
                    "area": [
                        latitude,
                        longitude,
                        latitude,
                        longitude,
                    ],
                    "format": "netcdf",
                },
                filepath,
            )


def convert_era5_data_wind(folder_path, output_csv="wind_all.csv"):
    all_dfs = []
    print(folder_path)
    # Find all NetCDF files in folder
    files = sorted(glob.glob(os.path.join(folder_path, "*.nc")))

    for file in files:
        print(f"Processing {file} ...")
        ds = xr.open_dataset(file, engine="h5netcdf")

        # Convert to pandas DataFrame
        df = ds[["u10", "v10"]].to_dataframe().reset_index()

        # Extract date components
        df["Year"] = df["valid_time"].dt.year
        df["Month"] = df["valid_time"].dt.month
        df["Day"] = df["valid_time"].dt.day
        df["Hour(utc)"] = df["valid_time"].dt.hour

        # Wind speed (m/s)
        df["301"] = np.sqrt(df["u10"] ** 2 + df["v10"] ** 2)

        # Wind direction (° from north, meteorological convention)
        df["365"] = (270 - np.degrees(np.arctan2(df["v10"], df["u10"]))) % 360

        # Keep only relevant columns
        df_final = df[["Year", "Month", "Day", "Hour(utc)", "301", "365"]]

        all_dfs.append(df_final)

        ds.close()

    # Merge all into one DataFrame
    df_all = pd.concat(all_dfs, ignore_index=True)

    # Save to CSV
    df_all.to_csv(output_csv, index=False, sep=";")
    print(f"Saved merged data to {output_csv}")

    return df_all


def download_era5_data_snow(
    latitude, longitude, year_start, year_end, output_folder="data/ERA5_snow"
):
    c = cdsapi.Client()

    # create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    years = [str(y) for y in range(year_start, year_end + 1)]
    months = [str(m).zfill(2) for m in range(1, 13)]  # all months
    days = [str(d).zfill(2) for d in range(1, 32)]  # all days
    times = ["00:00"]

    filename = f"snow_{year_start}_{year_end}.nc"
    filepath = os.path.join(output_folder, filename)
    print(f"Downloading {filename} for all years at once...")

    c.retrieve(
        "reanalysis-era5-single-levels",
        {
            "product_type": "reanalysis",
            "variable": ["skin_temperature"],
            # "variable": ["snow_albedo", "runoff", "skin_temperature"]
            "year": years,
            "month": months,
            "day": days,
            "time": times,
            "area": [latitude, longitude, latitude, longitude],
            "format": "netcdf",
        },
        filepath,
    )


def convert_era5_data_snow(folder_path, output_csv="snow_all.csv"):
    all_dfs = []
    print(folder_path)
    # Find all NetCDF files in folder
    files = sorted(glob.glob(os.path.join(folder_path, "*.nc")))

    for file in files:
        print(f"Processing {file} ...")
        ds = xr.open_dataset(file, engine="h5netcdf")

        print(ds)

        # Convert to pandas DataFrame
        df = ds[["skt"]].to_dataframe().reset_index()

        # Extract date components
        df["Year"] = df["valid_time"].dt.year
        df["Month"] = df["valid_time"].dt.month
        df["Day"] = df["valid_time"].dt.day
        df["Hour(utc)"] = df["valid_time"].dt.hour

        # Keep only relevant columns
        df_final = df[["Year", "Month", "Day", "Hour(utc)", "skt"]]

        all_dfs.append(df_final)

        print(df_final)

        ds.close()

    # Merge all into one DataFrame
    df_all = pd.concat(all_dfs, ignore_index=True)

    # Save to CSV
    df_all.to_csv(output_csv, index=False, sep=";")
    print(f"Saved merged data to {output_csv}")

    return df_all
