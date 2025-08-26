import pandas as pd
import matplotlib.pyplot as plt
from windrose import WindroseAxes
import numpy as np

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