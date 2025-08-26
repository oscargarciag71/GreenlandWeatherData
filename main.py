import pandas as pd
import matplotlib.pyplot as plt
import utils


# Read the CSV file
df = pd.read_csv("446500.csv", sep=";")
utils.plot_histograms(df)
