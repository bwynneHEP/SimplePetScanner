import pandas as pd
import matplotlib.pyplot as mpl

myColours=["grey", "grey", "red", "red", "gold", "gold", "yellowgreen", "yellowgreen", "deepskyblue", "deepskyblue", "navy", "navy", "hotpink", "hotpink"]

myLineStyles=['solid', 'dotted', 'solid', 'dotted', 'solid', 'dotted', 'solid', 'dotted', 'solid', 'dotted', 'solid', 'dotted', 'solid', 'dotted']

params = {'legend.fontsize': 15,
          'legend.title_fontsize': 15,
          'legend.loc': "upper left",
          'axes.labelsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15}
mpl.rcParams.update(params)

crystalsBlinded = [
    "BaF$_{2}$",
    "BGO",
    "CdWO$_4$",
    "CsF",
    "CsI",
    "CZT",
    "LSO",
    "LYSO",
    "NaI",
    "Perovskite 1",
    "Perovskite 2",
    "Perovskite 3",
	"Perovskite 4",
	"Perovskite 5", 
]

csvNames = ["maxNECR.csv", "trues.csv", "scatters.csv", "randoms.csv", "prompts.csv"]
labels = ["Max NECR [kcps]", "True rate [kcps]", "Scatter rate [kcps]", "Random rate [kcps]", "Prompt rate [kcps]"]
pdfNames = ["Siemens_maxNECR.pdf", "Siemens_trues.pdf", "Siemens_scatters.pdf", "Siemens_randoms.pdf", "Siemens_prompts.pdf"]

for ind, csvName in enumerate(csvNames):
    df = pd.read_csv(csvName)

    mpl.figure()
    mpl.gcf().set_size_inches(10, 10)
    index = range(1, 3, 1)
    x = df['0']
    for i in index :
        y = df[str(i)]*0.001
        mpl.plot(x, y, linewidth=4.0, color=myColours[i-1], label=crystalsBlinded[i-1], linestyle=myLineStyles[i-1])
    leglocation = 'upper right'
    if ind == 3 or ind == 4:
        leglocation = 'lower right'
    mpl.legend( crystalsBlinded, title="Crystals", loc=leglocation )
    mpl.xlabel( "Activity concentration [Bq/ml]" )
    mpl.ylabel( labels[ind] )
    mpl.title("Biograph Vision Quadra geometry", fontsize=20)
    mpl.savefig(pdfNames[ind])