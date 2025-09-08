import matplotlib.pyplot as mpl
import numpy as np
import copy

from CoincidenceGeneration import *
            
def TimelinesPlot( Timelines ):
    for i, decayTimes in enumerate( Timelines ):
        dummyY = []
        for time in decayTimes:
            dummyY.append( i )

        mpl.scatter( decayTimes, dummyY )

from SimulationDataset import *
import LegacyDatasetReader

RNG = np.random.default_rng(1)
print(RNG.random(10))
print(RNG.bit_generator.state)
N_CHANNELS = 20
BATCH_SIZE = 50
DECAY_RATES = [ float(i+1) for i in range( N_CHANNELS ) ]
decayTimes, timePeriod = TimeSeriesMultiChannel( BATCH_SIZE, DECAY_RATES, RNG )

# TimelinesPlot( decayTimes )

# mpl.xlabel( "Time [s]" )
# mpl.yticks( range( 0, N_CHANNELS ) )
# mpl.ylabel( "Decay channel" )
# mpl.gcf().set_size_inches( 10, 10 )
# mpl.show()

# just use any old existing file for this demo
dataSet = SimulationDataset( "hits.n1000000.SiemensBlock.1024mm.LinearF18.700mm.-y0mm.1234.csv", 1000000, RNG=RNG )
dataSetReader = LegacyDatasetReader.LegacyDatasetReader( dataSet )
inputData = [ dataSetReader ] * N_CHANNELS

# The "merged" part means merged across all channels, but for now show each separately
photonTimes = []
for i in range( N_CHANNELS ):
    photons = MergedPhotonStream( [decayTimes[i]], [inputData[i]], RNG )
    if len( photons ) > 0:
        photonTimes.append( photons[:, DATASET_TIME] )
    else:
        photonTimes.append( [] )

TimelinesPlot( photonTimes )

mpl.xlabel( "Time [ns]" )
mpl.yticks( range( 0, N_CHANNELS ) )
mpl.ylabel( "Decay channel" )
mpl.gcf().set_size_inches( 10, 10 )
mpl.show()

from matplotlib.patches import Rectangle

def CoincidenceBoxes( Coincidences, CoincidenceTimes, TimeWindow ):
    
    for i, coincidence in enumerate( Coincidences ):
        time = CoincidenceTimes[ i ]
        y = min( coincidence )
        height = max( coincidence ) - y
        y -= 0.1
        height += 0.2
        mpl.gca().add_patch( Rectangle( (time, y), TimeWindow, height, \
                                        linewidth=1,edgecolor='r',facecolor='none') )
# Not possible to reproduce previous plots exactly since RNG effects decays and photon sampling separately
RNG.bit_generator.state = np.random.default_rng(1).bit_generator.state
# print("RANDOM!!")
# print(RNG.random(10))
# print("STATE BEFORE")
# print(RNG.bit_generator.state)
coincidences = []
coincidenceTimes = []
photonTimes = []

inputDataCopy = copy.deepcopy(inputData)

for coincidence in GenerateCoincidences( BATCH_SIZE, DECAY_RATES, inputData, RNG, CoincidenceWindow=0.5E8, SimulationWindow=1E9, MultiWindow=False ):
    coincidenceTimes.append( coincidence[0, DATASET_TIME] )
    coincidences.append( [-0.5, 0.5] ) #photon channel info not preserved any more
    for photon in coincidence:
        photonTimes.append( photon[DATASET_TIME] )
photonTimes = [photonTimes]
# print("STATE AFTER")
# print(RNG.bit_generator.state)

TimelinesPlot( photonTimes )
CoincidenceBoxes( coincidences, coincidenceTimes, 0.5E8 )

mpl.xlabel( "Time [ns]" )
mpl.yticks( range( 0, 1 ) )
mpl.ylabel( "Decay channel" )
mpl.gcf().set_size_inches( 10, 1 )
# mpl.show()
mpl.savefig("normalcoincidences.pdf")

#Run once more to check reproducibility
print("#################################################")
print("#################################################")
print("#################################################")
RNG.bit_generator.state = np.random.default_rng(1).bit_generator.state
# print(RNG.random(10))
# print("STATE BEFORE")
# print(RNG.bit_generator.state)
coincidences = []
coincidenceTimes = []
photonTimes = []
for coincidence in GenerateDelayedCoincidences( BATCH_SIZE, DECAY_RATES, inputDataCopy, RNG, CoincidenceWindow=0.5E8, SimulationWindow=1E9, MultiWindow=False ):
    print("New delayed coincidence window")
    # coincidenceTimes.append( coincidence[0, DATASET_TIME] )
    # coincidences.append( [-0.5, 0.5] ) #photon channel info not preserved any more
    # for photon in coincidence:
    #     photonTimes.append( photon[DATASET_TIME] )
# photonTimes = [photonTimes]
# print("STATE AFTER")
# print(RNG.bit_generator.state)

# TimelinesPlot( photonTimes )
# CoincidenceBoxes( coincidences, coincidenceTimes, 0.5E8 )

# mpl.xlabel( "Time [ns]" )
# mpl.yticks( range( 0, 1 ) )
# mpl.ylabel( "Decay channel" )
# mpl.gcf().set_size_inches( 10, 1 )
# mpl.show()

# Not possible to reproduce previous plots exactly since RNG effects decays and photon sampling separately
# RNG.bit_generator.state = np.random.default_rng(1).bit_generator.state

# coincidences = []
# coincidenceTimes = []
# photonTimes = []
# for coincidence in GenerateCoincidences( BATCH_SIZE, DECAY_RATES, inputData, RNG, CoincidenceWindow=1E8, SimulationWindow=1E9, MultiWindow=True ):
#     coincidenceTimes.append( coincidence[0, DATASET_TIME] )
#     coincidences.append( [-0.5, 0.5] ) #photon channel info not preserved any more
#     for photon in coincidence:
#         photonTimes.append( photon[DATASET_TIME] )
# photonTimes = [photonTimes]

# TimelinesPlot( photonTimes )
# CoincidenceBoxes( coincidences, coincidenceTimes, 1E8 )

# mpl.xlabel( "Time [ns]" )
# mpl.yticks( range( 0, 1 ) )
# mpl.ylabel( "Decay channel" )
# mpl.gcf().set_size_inches( 10, 1 )
# mpl.show()