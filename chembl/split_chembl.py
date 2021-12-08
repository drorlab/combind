import sys
import pandas as pd

csv = sys.argv[1]
df = pd.read_csv(csv)

a = df['standard_value']
df[a < 10**2].to_csv(csv.replace('.csv', '_nM.csv'))
df[(a >= 10**2) & (a < 10**4)].to_csv(csv.replace('.csv', '_uM.csv'))
df[a >= 10**4].to_csv(csv.replace('.csv', '_mM.csv'))
