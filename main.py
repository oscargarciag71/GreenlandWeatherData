import pandas as pd
import matplotlib.pyplot as plt
import utils


# Read the CSV file
df = pd.read_csv("data/436000.csv", sep=";")
utils.plot_histograms(df)
utils.plot_wind_rose(df, 124)
