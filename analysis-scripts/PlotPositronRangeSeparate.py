import pandas as pd
import matplotlib.pyplot as mpl

colnames = ['event', 'd']
isotopeNames = ['$^{11}C$', '$^{13}N$', '$^{15}O$', '$^{18}F$', '$^{68}Ga$', '$^{82}Rb$']
isotopes = ['C11', 'N13', 'O15', 'F18', 'Ga68', 'Rb82']
maxRange = [3.8, 5.0, 8.0, 2.2, 9.0, 15.5]

def plotData():
    for i, isotope in enumerate(isotopes) :
        csvfileName = 'decays_' + isotope + '.csv'
        csvfile = pd.read_csv(csvfileName, delimiter=" ", names=colnames)

        mpl.rc('axes', axisbelow=True) #needed to put grid in the background
        ax = csvfile.plot.hist(column='d', range=[0, 16], bins=32, figsize=(6,6))
        ax.set_xlabel("Positron range [mm]")
        ax.get_legend().remove()

        #calculate mean positron range
        meanDist = csvfile['d'].mean()
        rounded = round(meanDist, 2)

        #calculate mean positron range for values that make sense... 
        indices = csvfile.index[csvfile['d']<maxRange[i]].tolist()
        meanDistCorr = csvfile['d'].iloc[indices].mean(axis=0)
        roundedCorrMean = round(meanDistCorr, 2)
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

if __name__ == '__main__':
    plotData()
