import subprocess

import numpy as np
import pandas as pd

# they will make a choice this program will run and feed it into main
filePath = input("Enter CSV File Path: ")

df = pd.read_csv(filePath)
df = df.dropna()
# now that we've dropped na we should be fine?

df.to_csv()

subprocess.run(["./build/OLS", "data/cleaned.csv"], check=True)
