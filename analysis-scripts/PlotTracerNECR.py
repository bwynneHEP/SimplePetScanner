import pandas as pd
import matplotlib.pyplot as mpl

params = {'legend.fontsize': 15,
          'legend.title_fontsize': 15,
          'legend.loc': "upper left",
          'axes.labelsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15}
mpl.rcParams.update(params)

dfC11 = pd.read_csv("tracerNECR_C11.csv")
dfN13 = pd.read_csv("tracerNECR_N13.csv")
dfO15 = pd.read_csv("tracerNECR_O15.csv")
dfF18 = pd.read_csv("tracerNECR_F18.csv")
dfGa68 = pd.read_csv("tracerNECR_Ga68.csv")
dfRb82 = pd.read_csv("tracerNECR_Rb82.csv")

dataframes = [dfC11, dfN13, dfO15, dfF18, dfGa68, dfRb82]
legItems = ['$^{11}$C', '$^{13}$N', '$^{15}$O', '$^{18}$F', '$^{68}$Ga', '$^{82}$Rb']
myColours=["grey", "red", "gold", "yellowgreen", "deepskyblue", "hotpink"]

mpl.figure()
mpl.gcf().set_size_inches(10, 10)

for i, df in enumerate(dataframes):
    x = df['activity']
    y = df['NECR']*0.001

    mpl.plot(x, y, linewidth=4.0, label=legItems[i], color=myColours[i])

mpl.legend( legItems, title="Isotopes", loc='upper right' )
mpl.xlabel( "Activity concentration [Bq/ml]" )
mpl.ylabel( "Max NECR [kcps]" )
mpl.title("Biograph Vision Quadra geometry with LSO crystals", fontsize=20)
mpl.savefig("tracerNECR_comparison_siemens_lso.pdf")