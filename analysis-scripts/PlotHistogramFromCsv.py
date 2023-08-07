import pandas as pd
import matplotlib.pyplot as mpl

csvfile = pd.read_csv('test.hits.n1000000.SiemensBlock.1024mm.LinearF18.700mm.1234.csv', delimiter=" ")

print(csvfile)
# mpl.clf()
csvfile.hist(column='R')
mpl.show()