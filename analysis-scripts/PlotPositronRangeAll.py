#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as mpl
import matplotlib
matplotlib.use('agg')

_numDataFiles = 6
_colours = (
    (231/255, 111/255,  81/255, 0.33),
    (244/255, 162/255,  97/255, 0.30),
    (233/255, 196/255, 106/255, 0.27),
    ( 42/255, 157/255, 143/255, 0.24),
    ( 38/255,  70/255,  83/255, 0.21),
    ( 74/255,  31/255,  85/255, 0.18),
)
_lineColours = (
    (231/255, 111/255,  81/255),
    (244/255, 162/255,  97/255),
    (233/255, 196/255, 106/255),
    ( 42/255, 157/255, 143/255),
    ( 38/255,  70/255,  83/255),
    ( 74/255,  31/255,  85/255),
)

colnames = ['d']
isotopeNames = ['$^{18}F$', '$^{11}C$', '$^{13}N$', '$^{15}O$', '$^{68}Ga$', '$^{82}Rb$']
isotopes = ['F18', 'C11', 'N13', 'O15', 'Ga68', 'Rb82']
maxRange = [2.2, 3.8, 5.0, 8.0, 9.0, 15.5]

def plotData():
    df = pd.DataFrame()

    for i, isotope in enumerate(isotopes) :
        print("Isotope = ", isotope)
        # The csv file has to be skimmed in order to not count the same positron twice
        # Easy way to do that: awk '!seen[$1]++' hits.csv > hits_skimmed.csv
        csvfileName = 'decays_' + isotope + '.csv'
        csvfile = pd.read_csv(csvfileName, delimiter=" ", names=colnames)

        #calculate the mean positron range
        indices = csvfile.index[csvfile['d']<maxRange[i]].tolist()
        meanDistCorr = csvfile['d'].iloc[indices].mean(axis=0)
        roundedCorrMean = round(meanDistCorr, 2)
        print("Mean of selected rows = ", roundedCorrMean)
        leg = isotopeNames[i] + ", <d> = " + str(roundedCorrMean)
        if i<4 :
            leg = isotopeNames[i] + ",   <d> = " + str(roundedCorrMean)
        df.insert(i, leg, csvfile['d'])

    ax1 = df.plot.hist(range=[0, 12], bins=32, figsize=(6, 6), color=_lineColours, histtype='step', lw=2)
    ax2 = df.plot.hist(range=[0, 12], bins=32, color=_colours, ax=ax1)
    
    ax1.set_xlabel("Positron range [mm]")

    handles, labels = mpl.gca().get_legend_handles_labels()
    mpl.legend(handles[0:df.columns.size], labels[0:df.columns.size])

    mpl.savefig('PositronRange-Siemens-LSO-All.pdf')

if __name__ == '__main__':
    plotData()
