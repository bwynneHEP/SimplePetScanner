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
    for time in range( 0, 700, 20 ):
        timeSec = float(time) * 60.0
        activity = TracerActivityAtTime( 1100E6, timeSec, isotope )

        activityList = []
        dataList = []
        # Add crystalActivity to the activityList only if the crystal material is radioactive
        if detectorMaterial == "LSO" or detectorMaterial == "LYSO" :
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

# First plot NECR vs. detectorLength for 4 phantom lengths for each isotope
# isotope = "Rb82"
detectorMaterial = "LSO"

def OneDetector( isotope, detectorLength, phantomLength, detectorMaterial, simulationWindow=1E-3, coincidenceWindow=4.7E-9 ):
    # Fix random seed for reproducibility, don't if you want variation
    # Has to be set in this method, not before, because this is where we enter the worker processes
    random.seed( detectorLength )

    tracer = "Linear" + isotope
    print("tracer = ", tracer)

    datasetSize = 100000

    tracerData = CreateDataset( detectorLength, "Siemens", phantomLength, tracer, datasetSize, 435.0, 585.0, detectorMaterial )
    crystalData = None
    crystalActivity = None
    #calculate crystal activity and create crystalDataset only if the crystal is radioactive
    if detectorMaterial == "LSO" or detectorMaterial == "LYSO" :
        crystalActivity= sqp.Lu176decaysInMass( sqp.DetectorMassLength( detectorLength, detectorMaterial ) )
        crystalData = CreateDataset( detectorLength, "Siemens", phantomLength, "Siemens", datasetSize, 435.0, 585.0, detectorMaterial )

    activityAtTimeSiemens, necrAtTimeSiemens, trueAtTimeSiemens, rPlusSAtTimeSiemens, scatterAtTimeSiemens, randomAtTimeSiemens = NECRatTime( isotope, tracerData, crystalData, crystalActivity, sqp.DetectorRadius(), phantomLength, detectorMaterial, simulationWindow, coincidenceWindow )
    return ( max( necrAtTimeSiemens ), sqp.DetectorDiscreteLength( detectorLength ) )

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

trialPhantoms = [ 300, 700, 1100, 1500 ]

isotopes = ["C11", "N13", "O15", "F18", "Ga68", "Rb82"]

for isotope in isotopes:
    detectorLengths = []
    maxNECRlines = []

    for phantomLength in trialPhantoms:
        detectorLengths, maxNECR = PeakNECRWithLengthMultiprocess( isotope, phantomLength, detectorMaterial, simulationWindow=1E-3, coincidenceWindow=4.7E-9 )
        maxNECRlines.append( maxNECR )

    for i, phantomLength in enumerate( trialPhantoms ):
        mpl.plot( detectorLengths, [val * 0.001 for val in maxNECRlines[i]], label=phantomLength, linewidth=4.0 )

    mpl.xlabel( "Detector length [mm]" )
    mpl.ylabel( "Max NECR [kcps]" )
    mpl.legend( trialPhantoms, title="Source length [mm]" )
    mpl.gcf().set_size_inches( 10, 10 )
    figName = 'NECRvsLength_' + isotope + '.pdf'
    mpl.savefig(figName)
    mpl.clf()

##########################################################
phantomLength = 700
datasetSize = 100000
crystalData = CreateDataset( 1024, "Siemens", phantomLength, "Siemens", datasetSize, siemensEmin, siemensEmax, detectorMaterial )
activityList = []
maxNECRlines = []
for isotope in isotopes:
    tracer = "Linear" + isotope
    print("tracer = ", tracer)
    tracerData = CreateDataset( 1024, "Siemens", phantomLength, tracer, datasetSize, siemensEmin, siemensEmax, detectorMaterial )

    activityAtTimeSiemens, necrAtTimeSiemens, trueAtTimeSiemens, rPlusSAtTimeSiemens, scatterAtTimeSiemens, randomAtTimeSiemens = NECRatTime( isotope, tracerData, crystalData, sqp.Lu176decaysInMass( sqp.DetectorMass(detectorMaterial) ), sqp.DetectorRadius(), phantomLength, detectorMaterial )
    activityList.append(activityAtTimeSiemens)
    maxNECRlines.append(necrAtTimeSiemens)
    promptAtTimeSiemens = [sum(n) for n in zip(trueAtTimeSiemens, scatterAtTimeSiemens, randomAtTimeSiemens)]

    mpl.plot( activityAtTimeSiemens, [val * 0.001 for val in necrAtTimeSiemens], linewidth=4.0, label="NECR" )
    mpl.plot( activityAtTimeSiemens, [val * 0.001 for val in promptAtTimeSiemens], linewidth=4.0, label="Prompts" )
    mpl.plot( activityAtTimeSiemens, [val * 0.001 for val in trueAtTimeSiemens], linewidth=4.0, label="Trues" )
    mpl.plot( activityAtTimeSiemens, [val * 0.001 for val in scatterAtTimeSiemens], linewidth=4.0, label="Scatter" )
    mpl.plot( activityAtTimeSiemens, [val * 0.001 for val in randomAtTimeSiemens], linewidth=4.0, label="Randoms" )
    legItems = ["NECR", "Prompts", "Trues", "Scatter", "Randoms"]
    mpl.legend(legItems)
    mpl.xlabel( "Activity [Bq/ml]" )
    mpl.ylabel( "Counts [kcps]" )
    mpl.gcf().set_size_inches(10, 10)
    figName = "TracerCounts_" + isotope + ".pdf"
    mpl.savefig(figName)
    mpl.clf()

mpl.clf()
for i, isotope in enumerate( isotopes ):
    mpl.plot( activityList[i], [val * 0.001 for val in maxNECRlines[i]], label=isotope, linewidth=4.0 )

mpl.legend( isotopes, title="Isotopes" )
mpl.xlim( [ 0, 20000 ] )
mpl.xlabel( "Activity [Bq/ml]" )
mpl.ylabel( "Counts [kcps]" )
mpl.gcf().set_size_inches(10, 10)
mpl.savefig("TracerNECR.pdf")

# ################################################

