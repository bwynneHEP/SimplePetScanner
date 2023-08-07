import sys
sys.path.append('../analysis/')
import SiemensQuadraProperties as sqp
from ActivityTools import *
from SimulationDataset import *

import multiprocessing as mp
mp.set_start_method('fork')
from multiprocessing import Pool
import random
import numpy as np
import matplotlib.pyplot as mpl
import sys

myColours=["royalblue", "darkorange", "yellowgreen", "hotpink", "orangered"]
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=myColours)

params = {'legend.fontsize': 15,
          'legend.title_fontsize': 15,
          'legend.loc': "lower right",
          'axes.labelsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15}
mpl.rcParams.update(params)

# Hanna: Set the simulation parameters globally
simulationWindow = 1E-3
coincidenceWindow = 4.7E-9
detectorLengths = [100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900]
phantomLength = 300
NDecays = 1000000

def NECRatTimeF18( tracerData, crystalData, crystalActivity, detectorRadius, phantomLength, detectorMaterial, simulationWindow=1E-3, coincidenceWindow=4.7E-9 ):
    # NEMU NU 2-2012 says 650mm window for 700mm phantom, so keep the same relation
    zWindow = (phantomLength - 50) / 2

    # get volume in cc
    phantomRadius = 20.3 / 2.0
    phantomVolume = phantomRadius * phantomRadius * math.pi * phantomLength / 10.0 # assume length in mm

    necrAtTime = []
    trueAtTime = []
    rPlusSAtTime = []
    activityAtTime = []
    for time in range( 0, 700, 20 ):
        timeSec = float(time) * 60.0
        activity = F18ActivityAtTime( 1100E6, timeSec )

        activityList = []
        dataList = []
        # Hanna: Add crystalActivity to the activityList only if the crystal material is radioactive
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
        activityAtTime.append( activity / phantomVolume )

    mpl.clf()
    return activityAtTime, necrAtTime, trueAtTime, rPlusSAtTime

def GetDataset(dirName, detectorLength, phantomLength, NDecays, detectorMaterial, energyMin=435.0, energyMax=585.0):

    random.seed( detectorLength )

    # Get tracer data
    tracerFileName = dirName + "hits.n" + str(NDecays) + ".SiemensBlock." + str(detectorLength) + "mm.LinearF18." + str(phantomLength) + "mm.1234.csv"
    tracerData = SimulationDataset(tracerFileName, NDecays, 435.0, 585.0)

    crystalActivity = None
    crystalData = None

    if detectorMaterial == "LSO" or detectorMaterial == "LYSO" :
        crystalActivity = sqp.Lu176decaysInMass( sqp.DetectorMassLength( detectorLength, detectorMaterial ) )
        # Get crystal data
        crystalFileName = dirName + "hits.n" + str(NDecays) + ".SiemensBlock." + str(detectorLength) + "mm.Siemens." + str(phantomLength) + "mm.1234.csv"
        crystalData = SimulationDataset(crystalFileName, NDecays, 435.0, 585.0)
    
    activityAtTimeSiemens, necrAtTimeSiemens, trueAtTimeSiemens, rPlusSAtTimeSiemens = NECRatTimeF18( tracerData, crystalData, crystalActivity, sqp.DetectorRadius(), phantomLength, detectorMaterial, simulationWindow, coincidenceWindow )
    return ( max( necrAtTimeSiemens ), sqp.DetectorDiscreteLength( detectorLength ) )

def PeakNECRWithLengthMultiprocess(dirName, phantomLength, detectorMaterial, simulationWindow=1E-3, coincidenceWindow=4.7E-9, processes=5 ):
    # Create the arguments for each process
    arguments = []
    for detectorLength in detectorLengths:
        arguments.append( (dirName, detectorLength, phantomLength, NDecays, detectorMaterial, simulationWindow, coincidenceWindow ) )

    # Launch a separate process for each detector length
    result = None
    with Pool( processes=processes ) as p:
        result = p.starmap( GetDataset, arguments )

    # Unpack the results
    realLengths = []
    maxNECR = []
    for entry in result:
        maxNECR.append( entry[0] )
        realLengths.append( entry[1] )

    # Useful when multiprocessing does not work, do not remove
    # for detectorLength in detectorLengths:
    #     result = GetDataset(dirName, detectorLength, phantomLength, NDecays, detectorMaterial, simulationWindow, coincidenceWindow)
    #     maxNECR.append(result[0])
    #     realLengths.append(result[1])
    
    return realLengths, maxNECR

mpl.figure()
legItems = ["LSO", "NaI", "BGO"]
items = ["LSO", "NaI", "BGO"]
dirNames = [
    "/Users/hannaborecka-bielska/Desktop/PET/pet-sim/sim/Data/LSO-1M/",
    "/Users/hannaborecka-bielska/Desktop/PET/pet-sim/sim/Data/NaI-1M/",
    "/Users/hannaborecka-bielska/Desktop/PET/pet-sim/sim/Data/BGO-1M/"
]
maxNECRlines = []
scannerTrueLengths = []

for i in range(len(legItems)):
    scannerTrueLengths, maxNECR = PeakNECRWithLengthMultiprocess(str(dirNames[i]), phantomLength, str(items[i]), simulationWindow=1E-3 )
    mpl.plot( scannerTrueLengths, maxNECR, linewidth=4.0, color= myColours[i], label= legItems[i] )

mpl.legend(title="Crystals" )
mpl.gcf().set_size_inches( 10, 10 )
mpl.xlabel( "Detector length [mm]" )
mpl.ylabel( "Max NECR [/sec]" )
pdfName = "compareNECR_crystal_" + str(phantomLength) + ".pdf"
mpl.savefig(pdfName)

    
    


