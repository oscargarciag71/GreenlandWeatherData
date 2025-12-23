import pandas as pd
import matplotlib.pyplot as plt
import utils


# Read the CSV file
df = pd.read_csv("data/446500.csv", sep=";")
utils.plot_histograms(df)

df = pd.read_csv("data/422100.csv", sep=";")
utils.plot_histograms(df)

df = pd.read_csv("data/447500.csv", sep=";")
utils.plot_histograms(df)

utils.plot_wind_rose(df, 304, 354, 254)
