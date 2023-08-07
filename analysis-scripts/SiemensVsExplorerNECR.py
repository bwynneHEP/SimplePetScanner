import sys
sys.path.append('../analysis/')
import SiemensQuadraProperties as sqp
import ExplorerProperties as ep
from ActivityTools import *
from SimulationDataset import *

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

# Calculate NECR versus time as an F18 linear source decays from an initial 1100MBq activity (plotted activity is divided by the cylindrical phantom volume in cc).

# Statistical fluctuation of the result can be reduced by increasing the simulationWindow parameter, at the cost of CPU time. The size of the simulated dataset can also be increased, but this has a lesser effect for counting experiments.

# The CreateDataset command will run the Geant4 simulation with the appropriate parameters for source and detector type and length. The requested number of events will be simulated, but only those within the specified energy window will be loaded for analysis. If a suitable dataset already exists it will be re-used without further simulation.

def NECRatTimeF18( tracerData, crystalData, crystalActivity, detectorRadius, phantomLength, simulationWindow=1E-2, coincidenceWindow=4.7E-9, zWindow=325.0 ):
    # get volume in cc
    phantomRadius = 20.3 / 2.0
    phantomVolume = phantomRadius * phantomRadius * math.pi * phantomLength / 10.0 # assume length in mm
    
    necrAtTime = []
    trueAtTime = []
    rPlusSAtTime = []
    activityAtTime = []
    scattersAtTime = []
    randomsAtTime = []

    for time in range( 0, 700, 20 ):
        timeSec = float(time) * 60.0
        activity = F18ActivityAtTime( 1100E6, timeSec )

        necr, true, rPlusS, Scatters, Randoms = DetectedCoincidences( [activity, crystalActivity], [tracerData, crystalData], simulationWindow, coincidenceWindow, detectorRadius, ZMin=-zWindow, ZMax=zWindow )
        necrAtTime.append( necr )
        trueAtTime.append( true )
        rPlusSAtTime.append( rPlusS )
        scattersAtTime.append( Scatters )
        randomsAtTime.append( Randoms )
        activityAtTime.append( activity / phantomVolume )
    
    return activityAtTime, necrAtTime, trueAtTime, rPlusSAtTime, scattersAtTime, randomsAtTime

# Simulation parameters
phantomLength = 700
datasetSize = 1000000
siemensEmin = 435.0
siemensEmax = 585.0
explorerEmin = 430.0
explorerEmax = 645.0

detectorMaterial = "LSO"
tracerData = CreateDataset( 1024, "Siemens", phantomLength, "LinearF18", datasetSize, siemensEmin, siemensEmax, detectorMaterial )
crystalData = CreateDataset( 1024, "Siemens", phantomLength, "Siemens", datasetSize, siemensEmin, siemensEmax, detectorMaterial )

activityAtTimeSiemens, necrAtTimeSiemens, trueAtTimeSiemens, rPlusSAtTimeSiemens, scatterAtTimeSiemens, randomAtTimeSiemens = NECRatTimeF18( tracerData, crystalData, sqp.Lu176decaysInMass( sqp.DetectorMass(detectorMaterial) ), sqp.DetectorRadius(), phantomLength )

mpl.clf()

detectorMaterial = "LYSO"
tracerData = CreateDataset( 1850, "Explorer", phantomLength, "LinearF18", datasetSize, explorerEmin, explorerEmax, detectorMaterial )
crystalData = CreateDataset( 1850, "Explorer", phantomLength, "Explorer", datasetSize, explorerEmin, explorerEmax, detectorMaterial )

activityAtTimeExplorer, necrAtTimeExplorer, trueAtTimeExplorer, rPlusSAtTimeExplorer, scatterAtTimeExplorer, randomAtTimeExplorer = NECRatTimeF18( tracerData, crystalData, ep.Lu176decaysInMass( ep.DetectorMass(detectorMaterial)), ep.DetectorRadius(), phantomLength )

mpl.clf()

labels = [ "Siemens NECR", "Explorer NECR", "Siemens True", "Explorer True", "Siemens R+S", "Explorer R+S" ]
mpl.plot( activityAtTimeSiemens, necrAtTimeSiemens, label=labels[0], linewidth=4.0 )
mpl.plot( activityAtTimeExplorer, necrAtTimeExplorer, label=labels[1], linewidth=4.0 )
mpl.plot( activityAtTimeSiemens, trueAtTimeSiemens, label=labels[2], linewidth=4.0 )
mpl.plot( activityAtTimeExplorer, trueAtTimeExplorer, label=labels[3], linewidth=4.0 )
mpl.plot( activityAtTimeSiemens, rPlusSAtTimeSiemens, label=labels[4], linewidth=4.0 )
mpl.plot( activityAtTimeExplorer, rPlusSAtTimeExplorer, label=labels[5], linewidth=4.0 )
mpl.legend( labels )
mpl.xlabel( "Activity [Bq/ml]" )
mpl.ylabel( "Counts [/sec]" )
mpl.gcf().set_size_inches(10, 10)
mpl.savefig("SiemensVsExplorer.pdf")
mpl.clf()

labels = [ "NECR", "True rate", "Scatter rate", "Random rate" ]
mpl.plot( activityAtTimeExplorer, necrAtTimeExplorer, label=labels[0], linewidth=4.0 )
mpl.plot( activityAtTimeExplorer, trueAtTimeExplorer, label=labels[1], linewidth=4.0 )
mpl.plot( activityAtTimeExplorer, scatterAtTimeExplorer, label=labels[2], linewidth=4.0 )
mpl.plot( activityAtTimeExplorer, randomAtTimeExplorer, label=labels[3], linewidth=4.0 )
mpl.legend( labels )
mpl.xlabel( "Activity [Bq/ml]" )
mpl.xlim( [ 0, 20000 ] )
mpl.ylim( [ 0, 8000000 ] )
mpl.ylabel( "Counts [/sec]" )
mpl.gcf().set_size_inches(10, 10)
mpl.savefig("Explorer_Counts.pdf")
mpl.clf()

mpl.plot( activityAtTimeSiemens, necrAtTimeSiemens, label=labels[0], linewidth=4.0 )
mpl.plot( activityAtTimeSiemens, trueAtTimeSiemens, label=labels[1], linewidth=4.0 )
mpl.plot( activityAtTimeSiemens, scatterAtTimeExplorer, label=labels[2], linewidth=4.0 )
mpl.plot( activityAtTimeSiemens, randomAtTimeExplorer, label=labels[3], linewidth=4.0 )
mpl.legend( labels )
mpl.xlabel( "Activity [Bq/ml]" )
mpl.xlim( [ 0, 40000 ] )
mpl.ylim( [ 0, 10000000 ] )
mpl.ylabel( "Counts [/sec]" )
mpl.gcf().set_size_inches(10, 10)
mpl.savefig("Siemens_Counts.pdf")
mpl.clf()