import pandas as pd
import matplotlib.pyplot as mpl

colnames = ['event', 'mod', 'energy', 'time', 'R', 'phi', 'z', 'd']
# isotopes = ['Rb82']
# isotopeNames = ['$^{82}Rb$']
# maxRange = [15.5]
isotopeNames = ['$^{11}C$', '$^{13}N$', '$^{15}O$', '$^{18}F$', '$^{68}Ga$', '$^{82}Rb$']
isotopes = ['C11', 'N13', 'O15', 'F18', 'Ga68', 'Rb82']
maxRange = [3.8, 5.0, 8.0, 2.2, 9.0, 15.5]

for i, isotope in enumerate(isotopes) :
    print("Isotope = ", isotope)
    # The csv file has to be skimmed in order to not count the same positron twice
    # Easy way to do that: awk '!seen[$1]++' hits.csv > hits_skimmed.csv
    csvfileName = 'hits_skimmed_Siemens_' + isotope + '.csv'
    csvfile = pd.read_csv(csvfileName, delimiter=" ", names=colnames)
    print(csvfile)

    mpl.rc('axes', axisbelow=True) #needed to put grid in the background
    ax = csvfile.plot.hist(column='d', range=[0, 16], bins=32, figsize=(6,6))
    ax.set_xlabel("Positron range [mm]")
    ax.get_legend().remove()

    #calculate mean positron range
    meanDist = csvfile['d'].mean()
    rounded = round(meanDist, 2)

    #calculate mean positron range for values that make sense... 
    indices = csvfile.index[csvfile['d']<maxRange[i]].tolist()
    # print("indices = ", indices)
    meanDistCorr = csvfile['d'].iloc[indices].mean(axis=0)
    outIndices = csvfile.index[csvfile['d']>maxRange[i]].tolist()
    outliers = csvfile['d'].iloc[outIndices]
    print("outliers = ", outliers)
    roundedCorrMean = round(meanDistCorr, 2)
    print("Mean of selected rows = ", roundedCorrMean)
    figtext = "Mean = " + str(roundedCorrMean)

    mpl.text(0.7, 0.9, figtext, fontsize='large', transform=ax.transAxes)
    mpl.grid(axis='both', color='0.95')
    title = 'Positron range for ' + isotopeNames[i]
    mpl.title(title)
    mpl.subplots_adjust(left=0.15)
    mpl.subplots_adjust(right=0.95)
    # mpl.show()
    pdfname = 'PositronRange-Siemens-LSO-' + isotope + '.pdf'
    mpl.savefig(pdfname)
    print('---------------------------------------')
