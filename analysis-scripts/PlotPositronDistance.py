import pandas as pd
import matplotlib.pyplot as mpl

colnames = ['event', 'mod', 'energy', 'time', 'R', 'phi', 'z', 'd']

# The csv file has to be skimmed in order to not count the same positron twice
# Easy way to do that: awk '!seen[$1]++' hits.csv > hits_skimmed.csv
csvfile = pd.read_csv('hits_skimmed.csv', delimiter=" ", names=colnames)
print(csvfile)

mpl.rc('axes', axisbelow=True) #needed to put grid in the background
ax = csvfile.plot.hist(column='d', range=[0, 16], bins=32, figsize=(6,6))
ax.set_xlabel("Positron distance [mm]")
ax.get_legend().remove()

#calculate mean positron range
meanDist = csvfile['d'].mean()
rounded = round(meanDist, 2)
figtext = "Mean = " + str(rounded)
print(figtext)

mpl.text(0.7, 0.9, figtext, fontsize='large', transform=ax.transAxes)
mpl.grid(axis='both', color='0.95')
mpl.title("Positron distance for $^{89}Zr$")
mpl.subplots_adjust(left=0.15)
mpl.subplots_adjust(right=0.95)
# mpl.show()
mpl.savefig("PositronRange-Siemens-LSO-Zr89.pdf")