import pandas as pd
import matplotlib.pyplot as plt


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
