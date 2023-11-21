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

# myColours=["royalblue", "darkorange", "yellowgreen", "hotpink", "orangered"]
myColours=["grey", "red", "gold", "yellowgreen", "deepskyblue", "hotpink"]
mpl.rcParams['axes.prop_cycle'] = mpl.cycler(color=myColours)

params = {'legend.fontsize': 15,
          'legend.title_fontsize': 15,
          'legend.loc': "lower right",
          'axes.labelsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15}
mpl.rcParams.update(params)

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

def GetDataset(dirName, detectorLength, phantomLength, NDecays, detectorMaterial, Emin, Emax, simulationWindow, coincidenceWindow):

    random.seed( detectorLength )

    # Get tracer data
    tracerFileName = dirName + "hits.n" + str(NDecays) + ".SiemensBlock." + str(detectorLength) + "mm.LinearF18." + str(phantomLength) + "mm.1234.csv"
    tracerData = SimulationDataset(tracerFileName, NDecays, Emin, Emax)

    crystalActivity = None
    crystalData = None

    if detectorMaterial == "LSO" or detectorMaterial == "LYSO" :
        crystalActivity = sqp.Lu176decaysInMass( sqp.DetectorMassLength( detectorLength, detectorMaterial ) )
        # Get crystal data
        crystalFileName = dirName + "hits.n" + str(NDecays) + ".SiemensBlock." + str(detectorLength) + "mm.Siemens." + str(phantomLength) + "mm.1234.csv"
        crystalData = SimulationDataset(crystalFileName, NDecays, Emin, Emax)
    
    activityAtTimeSiemens, necrAtTimeSiemens, trueAtTimeSiemens, rPlusSAtTimeSiemens = NECRatTimeF18( tracerData, crystalData, crystalActivity, sqp.DetectorRadius(), phantomLength, detectorMaterial, simulationWindow, coincidenceWindow )
    return ( max( necrAtTimeSiemens ), sqp.DetectorDiscreteLength( detectorLength ) )

def PeakNECRWithLengthMultiprocess(dirName, phantomLength, NDecays, detectorMaterial, Emin, Emax, simulationWindow, coincidenceWindow, processes=10 ):
    # Create the arguments for each process
    arguments = []
    for detectorLength in detectorLengths:
        arguments.append( (dirName, detectorLength, phantomLength, NDecays, detectorMaterial, Emin, Emax, simulationWindow, coincidenceWindow ) )

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
legItems = ["BGO", "CdWO$_{4}$", "LSO", "Perovskite 1", "Perovskite 3", "Perovskite 5"]
items = ["BGO", "CdWO4", "LSO", "CsPbBr3", "FAPbI3", "MAPbI3"]
dirNames = [
    "/Users/hannaborecka-bielska/Desktop/PET/pet-sim/sim/SimplePetScanner/analysis-scripts-bgo/",
    "/Users/hannaborecka-bielska/Desktop/PET/pet-sim/sim/SimplePetScanner/analysis-scripts-cdwo4/",
    "/Users/hannaborecka-bielska/Desktop/PET/pet-sim/sim/SimplePetScanner/analysis-scripts-lso/",
    "/Users/hannaborecka-bielska/Desktop/PET/pet-sim/sim/SimplePetScanner/analysis-scripts-cspbbr3/",
    "/Users/hannaborecka-bielska/Desktop/PET/pet-sim/sim/SimplePetScanner/analysis-scripts-fapbi3/",
    "/Users/hannaborecka-bielska/Desktop/PET/pet-sim/sim/SimplePetScanner/analysis-scripts-mapbi3/",
]
# Hanna: Set the simulation parameters globally
simulationWindow = 1E-3
coincidenceWindow = 4.7E-9
detectorLengths = [100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900]
phantomLength = 700
NDecays = 1000000
Emin = 435.0
Emax = 585.0

maxNECRlines = []
scannerTrueLengths = []

for i in range(len(legItems)):
    scannerTrueLengths, maxNECR = PeakNECRWithLengthMultiprocess(str(dirNames[i]), phantomLength, NDecays, str(items[i]), Emin, Emax, simulationWindow, coincidenceWindow )
    mpl.plot( scannerTrueLengths, [val * 0.001 for val in maxNECR], linewidth=4.0, color= myColours[i], label= legItems[i] )

mpl.legend(title="Crystals", loc='upper left' )
mpl.gcf().set_size_inches( 10, 10 )
mpl.xlabel( "Detector length [mm]" )
mpl.ylabel( "Max NECR [kcps]" )
mpl.title("Biograph Vision Quadra geometry", fontsize=20)
pdfName = "compareNECR_crystal_" + str(phantomLength) + ".pdf"
mpl.savefig(pdfName)

    
    


