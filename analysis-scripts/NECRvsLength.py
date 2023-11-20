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
# myColours=[royalblue, darkorange, yellowgreen, hotpink, orangered]
myColours = []
for colour in mpl.cm.viridis( np.linspace( 0.05, 0.8, 4 ) ):
    myColours.append( colour )
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

# Investigate ideal detector length with Siemens geometry
def OneDetector( detectorLength, phantomLength, detectorMaterial, nevents, Emin, Emax,  simulationWindow=1E-3, coincidenceWindow=4.7E-9 ):
    # Fix random seed for reproducibility, don't if you want variation
    # Has to be set in this method, not before, because this is where we enter the worker processes
    random.seed( detectorLength )

    tracerData = CreateDataset( detectorLength, "Siemens", phantomLength, "LinearF18", nevents, Emin, Emax, detectorMaterial )
    crystalData = None
    crystalActivity = None
    #calculate crystal activity and create crystalDataset only if the crystal is radioactive
    if detectorMaterial == "LSO" or detectorMaterial == "LYSO" :
        crystalActivity= sqp.Lu176decaysInMass( sqp.DetectorMassLength( detectorLength, detectorMaterial ) )
        crystalData = CreateDataset( detectorLength, "Siemens", phantomLength, "Siemens", nevents, Emin, Emax, detectorMaterial )

    activityAtTimeSiemens, necrAtTimeSiemens, trueAtTimeSiemens, rPlusSAtTimeSiemens = NECRatTimeF18( tracerData, crystalData, crystalActivity, sqp.DetectorRadius(), phantomLength, detectorMaterial, simulationWindow, coincidenceWindow )
    return ( max( necrAtTimeSiemens ), sqp.DetectorDiscreteLength( detectorLength ) )

def PeakNECRWithLengthMultiprocess( phantomLength, detectorMaterial, nevents, Emin, Emax, simulationWindow=1E-3, coincidenceWindow=4.7E-9, processes=10 ):
    # Create the arguments for each process
    detectorLengths = [ 100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900 ]
    arguments = []
    realLengths = []
    maxNECR = []
    for detectorLength in detectorLengths:
        arguments.append( ( detectorLength, phantomLength, detectorMaterial, nevents, Emin, Emax, simulationWindow, coincidenceWindow ) )
        # result = OneDetector(detectorLength, phantomLength, detectorMaterial, simulationWindow, coincidenceWindow)
        # maxNECR.append(result[0])
        # realLengths.append(result[1])
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

detectorLengths = []
maxNECRlines = []
trialPhantoms = [ 300, 700, 1100, 1500 ]

# Simulation Settings -- check carefully!
detectorMaterial = "LSO"
nevents = 1000000
Emin = 435.0
Emax = 585.0
simulationWindow = 1E-3
coincidenceWindow = 4.7E-9

for phantomLength in trialPhantoms:
    detectorLengths, maxNECR = PeakNECRWithLengthMultiprocess( phantomLength, detectorMaterial, nevents, Emin, Emax, simulationWindow, coincidenceWindow )
    maxNECRlines.append( maxNECR )

for i, phantomLength in enumerate( trialPhantoms ):
    mpl.plot( detectorLengths, [val * 0.001 for val in maxNECRlines[i]], label=phantomLength, linewidth=4.0 )

mpl.xlabel( "Detector length [mm]" )
mpl.ylabel( "Max NECR [kcps]" )
mpl.legend( trialPhantoms, title="Source length [mm]" )
mpl.gcf().set_size_inches( 10, 10 )
mpl.savefig('NECRvsLength.pdf')

def scale( y0, y1, x1 ):
    m = (y1 - y0) / (x1 - detectorLengths[0])
    y2 = m * ( detectorLengths[-1] - x1 )
    y2 += y1
    return [ y0, y2 ]

def plotThreeParts( xs, ys, i, label, colour ):
    theplot, = mpl.plot( xs, ys, label=label, color=colour, linewidth=4.0 )
    mpl.plot( [xs[0], xs[-1]], scale(ys[0], ys[i], xs[i]), label="_", color=colour, linestyle="dashed" )
    mpl.scatter( [xs[i]], [ys[i]], label="_", color=colour )
    return theplot

legendEntries = []
for i, trialPhantom in enumerate( trialPhantoms ):

    # find the detector length index closest to the source length
    detectorLengthIndex = 0
    while detectorLengths[ detectorLengthIndex ] < trialPhantom:
        detectorLengthIndex += 1

    theplot = plotThreeParts( detectorLengths, [val * 0.001 for val in maxNECRlines[i]], detectorLengthIndex, trialPhantoms[i], myColours[i] )
    legendEntries.append( theplot )

mpl.xlabel( "Detector length [mm]")
mpl.ylabel( "Max NECR [kcps]")
mpl.legend( legendEntries, trialPhantoms, title="Source length [mm]" )
mpl.gcf().set_size_inches(10,10)
mpl.savefig('NECRvsLength_wLine.pdf')
