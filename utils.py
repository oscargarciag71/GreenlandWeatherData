import pandas as pd
import matplotlib.pyplot as plt
from windrose import WindroseAxes
import numpy as np


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
        subset["101"].dropna().hist(bins=30, ax=axes[0], edgecolor="black")
        axes[0].set_title("Mean air temperature")
        axes[0].set_xlabel("ºC")
        # axes[0].set_ylabel("Frequency")

        # Histogram for column 301
        subset["301"].dropna().hist(bins=30, ax=axes[1], edgecolor="black")
        axes[1].set_title("Mean wind speed")
        axes[1].set_xlabel("m/s")
        # axes[1].set_ylabel("Frequency")

        plt.tight_layout(rect=[0, 0, 1, 0.95])
        plt.show()

def plot_wind_rose(df, marker_angle, kite_angle_1, kite_angle_2):
    # Months of interest
    months = {5: "May", 6: "June", 7: "July"}

    fig = plt.figure(figsize=(18, 6))

    # Create 3 windrose axes manually
    for i, (month, name) in enumerate(months.items()):
        df_month = df[df['Month'] == month][['301','365']].dropna()
        ws = df_month['301'].values
        wd = df_month['365'].values

        # Define axes rectangle: left, bottom, width, height
        rect = [0.05 + i*0.32, 0.1, 0.3, 0.8]  # adjust spacing between subplots
        ax = WindroseAxes(fig, rect)
        fig.add_axes(ax)
        
        ax.bar(wd, ws, opening=0.8, edgecolor="white", bins=[0,2,4,6,8,10], normed=True)
        ax.set_title(f"Average wind speed and direction for {name}", fontsize=14)
        
        # Add radial marker at 124°
        theta = np.deg2rad(marker_angle)
        ax.plot([theta, theta], [0, ax.get_rmax()], "r--", lw=2)

        theta_1 = np.deg2rad(kite_angle_1)
        ax.plot([theta_1, theta_1], [0, ax.get_rmax()], "b--", lw=2)

        theta_2 = np.deg2rad(kite_angle_2)
        ax.plot([theta_2, theta_2], [0, ax.get_rmax()], "b--", lw=2)
        
    # Add a single legend for all subplots
    ax.set_legend(title="Wind speed (m/s)")

    plt.show()
