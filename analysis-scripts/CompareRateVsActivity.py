import sys
sys.path.append('../analysis/')
import SiemensQuadraProperties as sqp
import ExplorerProperties as ep
from ActivityTools import *
from SimulationDataset import *

import multiprocessing as mp
mp.set_start_method('fork')
from multiprocessing import Pool
import matplotlib.pyplot as mpl
import pandas as pd

myColours=["grey", "grey", "red", "red", "gold", "gold", "yellowgreen", "yellowgreen", "deepskyblue", "deepskyblue", "navy", "navy", "hotpink", "hotpink"]
myLineStyles=['solid', 'dotted', 'solid', 'dotted', 'solid', 'dotted', 'solid', 'dotted', 'solid', 'dotted', 'solid', 'dotted', 'solid', 'dotted']
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=myColours)

params = {'legend.fontsize': 15,
          'legend.title_fontsize': 15,
          'legend.loc': "upper left",
          'axes.labelsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15}
mpl.rcParams.update(params)

#Fix random seed for reproducibility, or skip this to allow the results to vary
import random
random.seed( 1234 )

def NECRatTime( isotope, tracerData, crystalData, crystalActivity, detectorRadius, phantomLength, detectorMaterial, simulationWindow=1E-3, coincidenceWindow=4.7E-9 ):
    # NEMU NU 2-2012 says 650mm window for 700mm phantom, so keep the same relation
    zWindow = (phantomLength - 50) / 2

    # get volume in cc
    phantomRadius = 20.3 / 2.0
    phantomVolume = phantomRadius * phantomRadius * math.pi * phantomLength / 10.0 # assume length in mm

    necrAtTime = []
    trueAtTime = []
    rPlusSAtTime = []
    scatterAtTime = []
    randomAtTime = []
    activityAtTime = []
    timeRange = range( 0, 700, 20 )
    if isotope == "O15" or isotope == "Rb82":
        timeRange = range( 0, 30, 1 )
    if isotope == "C11":
        timeRange = range(0, 300, 5)
    if isotope == "N13":
        timeRange = range(0, 150, 2)
    # for time in range( 0, 700, 20 ):
    for time in timeRange:
        timeSec = float(time) * 60.0
        activity = TracerActivityAtTime( 1100E6, timeSec, isotope )

        activityList = []
        dataList = []
        # Add crystalActivity to the activityList only if the crystal material is radioactive
        if detectorMaterial == "LSO" or detectorMaterial == "LYSO" or detectorMaterial == "eCsI":
            activityList = [activity, crystalActivity]
            dataList = [tracerData, crystalData]
        else :
            activityList = [activity]
            dataList = [tracerData]
        necr, true, rPlusS, scatters, randoms = DetectedCoincidences( activityList, dataList, simulationWindow, coincidenceWindow, detectorRadius, ZMin=-zWindow, ZMax=zWindow )
        necrAtTime.append( necr )
        trueAtTime.append( true )
        rPlusSAtTime.append( rPlusS )
        scatterAtTime.append( scatters )
        randomAtTime.append( randoms )
        activityAtTime.append( activity / phantomVolume )

    mpl.clf()
    return activityAtTime, necrAtTime, trueAtTime, rPlusSAtTime, scatterAtTime, randomAtTime

# Simulation parameters
# phantomLength = 700
siemensEmin = 435.0
siemensEmax = 585.0
explorerEmin = 430.0
explorerEmax = 645.0

def OneDetector( isotope, detectorLength, phantomLength, detectorMaterial, simulationWindow=1E-3, coincidenceWindow=4.7E-9 ):
    # Fix random seed for reproducibility, don't if you want variation
    # Has to be set in this method, not before, because this is where we enter the worker processes
    random.seed( detectorLength )

    tracer = "Linear" + isotope
    print("tracer = ", tracer)

    datasetSize = 1000000

    tracerData = CreateDataset( detectorLength, "Siemens", phantomLength, tracer, datasetSize, siemensEmin, siemensEmax, detectorMaterial )
    crystalData = None
    crystalActivity = None
    #calculate crystal activity and create crystalDataset only if the crystal is radioactive
    if detectorMaterial == "LSO" or detectorMaterial == "LYSO" or detectorMaterial == "eCsI":
        if detectorMaterial == "eCsI":
            crystalActivity= sqp.Cs137decaysInMass( sqp.DetectorMassLength( detectorLength, detectorMaterial ) )
        else:
            crystalActivity= sqp.Lu176decaysInMass( sqp.DetectorMassLength( detectorLength, detectorMaterial ) )
        crystalData = CreateDataset( detectorLength, "Siemens", phantomLength, "Siemens", datasetSize, siemensEmin, siemensEmax, detectorMaterial )
    

    activityAtTimeSiemens, necrAtTimeSiemens, trueAtTimeSiemens, rPlusSAtTimeSiemens, scatterAtTimeSiemens, randomAtTimeSiemens = NECRatTime( isotope, tracerData, crystalData, crystalActivity, sqp.DetectorRadius(), phantomLength, detectorMaterial, simulationWindow, coincidenceWindow )
    return ( max( necrAtTimeSiemens ), sqp.DetectorDiscreteLength( detectorLength ) )

def GetDetDataset(dirName, isotope, detectorLength, phantomLength, NDecays, detectorMaterial, energyMin=siemensEmin, energyMax=siemensEmax, simulationWindow=1E-3, coincidenceWindow=4.7E-9):

    random.seed( detectorLength )

    crystalActivity = None
    crystalData = None

    if detectorMaterial == "LSO" or detectorMaterial == "LYSO" :
        crystalActivity = sqp.Lu176decaysInMass( sqp.DetectorMassLength( detectorLength, detectorMaterial ) )
        # Get crystal data
        crystalFileName = dirName + "hits.n" + str(NDecays) + ".SiemensBlock." + str(detectorLength) + "mm.Siemens." + str(phantomLength) + "mm.1234.csv"
        crystalData = SimulationDataset(crystalFileName, NDecays, siemensEmin, siemensEmax)
    
    return crystalData

def GetTracerDataset(dirName, isotope, detectorLength, phantomLength, NDecays, detectorMaterial, energyMin=siemensEmin, energyMax=siemensEmax, simulationWindow=1E-3, coincidenceWindow=4.7E-9):

    random.seed( detectorLength )
    tracerFileName = dirName + "hits.n" + str(NDecays) + ".SiemensBlock." + str(detectorLength) + "mm.LinearF18." + str(phantomLength) + "mm.1234.csv"
    print("tracerFileName = ", tracerFileName)
    tracerData = SimulationDataset(tracerFileName, NDecays, siemensEmin, siemensEmax)
    print("tracerData = ", tracerData)
    return tracerData


def PeakNECRWithLengthMultiprocess( isotope, phantomLength, detectorMaterial, simulationWindow=1E-3, coincidenceWindow=4.7E-9, processes=10 ):
    # Create the arguments for each process
    detectorLengths = [ 100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900 ]
    arguments = []
    realLengths = []
    maxNECR = []
    for detectorLength in detectorLengths:
        arguments.append( ( isotope, detectorLength, phantomLength, detectorMaterial, simulationWindow, coincidenceWindow ) )
    
    # Launch a separate process for each detector length
    result = None
    with Pool( processes=processes ) as p:
        result = p.starmap( OneDetector, arguments )

    # Unpack the results
    realLengths = []
    maxNECR = []
    for entry in result:
        maxNECR.append( entry[0] )
        realLengths.append( entry[1] )
    return realLengths, maxNECR

isotope = "F18"
simWin = 1E-2
phantomLength = 700
datasetSize = 1000000
crystalData = None
crystalActivity = None
detectorLength = 1024

#alphabetical order
crystals = [
    "BaF2",
    "BGO",
    "CdWO4",
    "CsF",
    "CsI",
    "CZT",
    "LSO",
    "LYSO",
    "NaI",
    "CsPbBr3",
    "Cs2AgBiBr6",
    "FAPbI3",
    "MAPbBr3",
	"MAPbI3",
]

crystalsBlinded = [
    "BaF2",
    "BGO",
    "CdWO4",
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

activityList = []
maxNECRlines = []
truelines = []
scatterlines = []
randomlines = []
promptlines = []

for detectorMaterial in crystals:

    #calculate crystal activity and create crystalDataset only if the crystal is radioactive
    #if crystal is Lutetium-based
    if detectorMaterial == "LSO":
        crystalActivity= sqp.Lu176decaysInMass( sqp.DetectorMassLength( detectorLength, detectorMaterial ) )
        crystalData = GetDetDataset("/Users/hannaborecka-bielska/Desktop/PET/pet-sim/sim/SimplePetScanner/analysis-scripts-lso/", isotope, detectorLength, phantomLength, datasetSize, detectorMaterial, siemensEmin, siemensEmax, simWin)
    if detectorMaterial == "LYSO":
        crystalActivity= sqp.Lu176decaysInMass( sqp.DetectorMassLength( detectorLength, detectorMaterial ) )
        crystalData = GetDetDataset("/Users/hannaborecka-bielska/Desktop/PET/pet-sim/sim/SimplePetScanner/analysis-scripts-lyso/", isotope, detectorLength, phantomLength, datasetSize, detectorMaterial, siemensEmin, siemensEmax, simWin)

    tracer = "LinearF18"
    tracerDirName = "/Users/hannaborecka-bielska/Desktop/PET/pet-sim/sim/SimplePetScanner/analysis-scripts-" + detectorMaterial.lower() + "/"
    tracerData = GetTracerDataset(tracerDirName, isotope, detectorLength, phantomLength, datasetSize, detectorMaterial, energyMin=siemensEmin, energyMax=siemensEmax, simulationWindow=1E-3, coincidenceWindow=4.7E-9)

    activityAtTimeSiemens, necrAtTimeSiemens, trueAtTimeSiemens, rPlusSAtTimeSiemens, scatterAtTimeSiemens, randomAtTimeSiemens = NECRatTime( "F18", tracerData, crystalData, crystalActivity, sqp.DetectorRadius(), phantomLength, detectorMaterial, simulationWindow=simWin )
    promptAtTimeSiemens = [sum(n) for n in zip(trueAtTimeSiemens, scatterAtTimeSiemens, randomAtTimeSiemens)]
    activityList.append(activityAtTimeSiemens)
    maxNECRlines.append(necrAtTimeSiemens)
    truelines.append(trueAtTimeSiemens)
    scatterlines.append(scatterAtTimeSiemens)
    randomlines.append(randomAtTimeSiemens)
    promptlines.append(promptAtTimeSiemens)

mpl.clf()
for i, detectorMaterial in enumerate( crystalsBlinded ):
    mpl.plot( activityList[i], [val * 0.001 for val in maxNECRlines[i]], label=detectorMaterial, linewidth=4.0, color=myColours[i], linestyle=myLineStyles[i] )

dfNECR = pd.DataFrame(list(zip(*[activityList[0], maxNECRlines[0], maxNECRlines[1]]
                                 , maxNECRlines[2], maxNECRlines[3], maxNECRlines[4], maxNECRlines[5], maxNECRlines[6], maxNECRlines[7], maxNECRlines[8], maxNECRlines[9], maxNECRlines[10], maxNECRlines[11], maxNECRlines[12], maxNECRlines[13]] )))
print("dfNECR = ", dfNECR)
dfNECR.to_csv('maxNECR.csv', index=False)

mpl.legend( crystalsBlinded, title="Crystals", loc='upper right' )
mpl.xlabel( "Activity [Bq/ml]" )
mpl.ylabel( "Max NECR [kcps]" )
mpl.gcf().set_size_inches(10, 10)
mpl.savefig("CompareRateVsActivity_NECR.pdf")

mpl.clf()
for i, detectorMaterial in enumerate( crystalsBlinded ):
    mpl.plot( activityList[i], [val * 0.001 for val in truelines[i]], label=detectorMaterial, linewidth=4.0, color=myColours[i], linestyle=myLineStyles[i] )

dfTrues = pd.DataFrame(list(zip(*[activityList[0], truelines[0], truelines[1]]
                                  , truelines[2], truelines[3], truelines[4], truelines[5], truelines[6], truelines[7], truelines[8], truelines[9], truelines[10], truelines[11], truelines[12], truelines[13]] )))
print("dfTrues = ", dfTrues)
dfTrues.to_csv('trues.csv', index=False)

mpl.legend( crystalsBlinded, title="Crystals", loc='upper right' )
mpl.xlabel( "Activity [Bq/ml]" )
mpl.ylabel( "True rate [kcps]" )
mpl.gcf().set_size_inches(10, 10)
mpl.savefig("CompareRateVsActivity_trues.pdf")

mpl.clf()
for i, detectorMaterial in enumerate( crystalsBlinded ):
    mpl.plot( activityList[i], [val * 0.001 for val in scatterlines[i]], label=detectorMaterial, linewidth=4.0, color=myColours[i], linestyle=myLineStyles[i] )

dfScatters = pd.DataFrame(list(zip(*[activityList[0], scatterlines[0], scatterlines[1]]
                                     , scatterlines[2], scatterlines[3], scatterlines[4], scatterlines[5], scatterlines[6], scatterlines[7], scatterlines[8], scatterlines[9], scatterlines[10], scatterlines[11], scatterlines[12], scatterlines[13]] )))
print("dfScatters = ", dfScatters)
dfScatters.to_csv('scatters.csv', index=False)

mpl.legend( crystalsBlinded, title="Crystals", loc='upper right' )
mpl.xlabel( "Activity [Bq/ml]" )
mpl.ylabel( "Scatter rate [kcps]" )
mpl.gcf().set_size_inches(10, 10)
mpl.savefig("CompareRateVsActivity_scatters.pdf")

mpl.clf()
for i, detectorMaterial in enumerate( crystalsBlinded ):
    mpl.plot( activityList[i], [val * 0.001 for val in randomlines[i]], label=detectorMaterial, linewidth=4.0, color=myColours[i], linestyle=myLineStyles[i] )

dfRandoms = pd.DataFrame(list(zip(*[activityList[0], randomlines[0], randomlines[1]]
                                    , randomlines[2], randomlines[3], randomlines[4], randomlines[5], randomlines[6], randomlines[7], randomlines[8], randomlines[9], randomlines[10], randomlines[11], randomlines[12], randomlines[13]] )))
print("dfRandoms = ", dfRandoms)
dfRandoms.to_csv('randoms.csv', index=False)

mpl.legend( crystalsBlinded, title="Crystals", loc='lower right' )
mpl.xlabel( "Activity [Bq/ml]" )
mpl.ylabel( "Random rate [kcps]" )
mpl.gcf().set_size_inches(10, 10)
mpl.savefig("CompareRateVsActivity_randoms.pdf")

mpl.clf()
for i, detectorMaterial in enumerate( crystalsBlinded ):
    mpl.plot( activityList[i], [val * 0.001 for val in promptlines[i]], label=detectorMaterial, linewidth=4.0, color=myColours[i], linestyle=myLineStyles[i] )

dfPrompts = pd.DataFrame(list(zip(*[activityList[0], promptlines[0], promptlines[1]]
                                    , promptlines[2], promptlines[3], promptlines[4], promptlines[5], promptlines[6], promptlines[7], promptlines[8], promptlines[9], promptlines[10], promptlines[11], promptlines[12], promptlines[13]] )))
print("dfPrompts = ", dfPrompts)
dfPrompts.to_csv('prompts.csv', index=False)

mpl.legend( crystalsBlinded, title="Crystals", loc='lower right' )
mpl.xlabel( "Activity [Bq/ml]" )
mpl.ylabel( "Prompt rate [kcps]" )
mpl.gcf().set_size_inches(10, 10)
mpl.savefig("CompareRateVsActivity_prompts.pdf")

print("maxNECRlines = ", maxNECRlines)
print("")
print("truelines = ", truelines)
print("")
print("scatterlines = ", scatterlines)
print("")
print("randomlines = ", randomlines)
print("")
print("promptlines = ", promptlines)