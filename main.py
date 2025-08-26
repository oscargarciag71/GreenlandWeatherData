import pandas as pd
import matplotlib.pyplot as plt

# Read the CSV file
df = pd.read_csv("446500.csv", sep=";")

# Plot histogram of column 101
plt.figure(figsize=(8, 5))
df["101"].hist(bins=30, edgecolor="black")
plt.title("Histogram of Column 101")
plt.xlabel("Value")
plt.ylabel("Frequency")
plt.grid(False)
plt.show()
