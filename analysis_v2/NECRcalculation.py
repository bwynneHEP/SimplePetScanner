import matplotlib.pyplot as plt
import SimulationDataset as sd


# Similar to NEMA standard but uses a 1D histogram
# Assumes source is at z-axis
def NECRFromHistogram( bins, values, SimulationWindow ):

    # Assume one bin edge is on the 20mm boundary
    lowValue = 0.0
    lowBin = 0
    highValue = 0.0
    highBin = 0
    histogramTotal = 0.0
    centralTotal = 0.0
    for binIndex in range( len( bins ) - 1 ):
        if bins[binIndex] == -20.0:
            lowValue = ( values[binIndex] + values[binIndex-1] ) / 2.0
            lowBin = binIndex
        elif bins[binIndex] == 20.0:
            highValue = ( values[binIndex] + values[binIndex-1] ) / 2.0
            highBin = binIndex
    
        if bins[binIndex] >= -120.0 and bins[binIndex] < 120.0:
            histogramTotal += values[binIndex]
        if bins[binIndex] >= -20.0 and bins[binIndex] < 20.0:
            centralTotal += values[binIndex]
            
    centralBackground = ( lowValue + highValue ) / 2.0
    centralBackground *= ( highBin - lowBin ) # number of bins in the region
    
    trues = centralTotal - centralBackground
    rPlusS = histogramTotal - trues
    necr = 0.0
    if histogramTotal > 0.0:
        necr = trues * trues / histogramTotal

    return necr/SimulationWindow, trues/SimulationWindow, rPlusS/SimulationWindow


# Run a coincidence generator and calculate NECR
def NECRfromGenerator( Generator, DetectorRadius, SimulationWindow, PairMode="Exclusive", ZMin=0.0, ZMax=0.0 ):

    hitRadii = []
    hitRadiiTruth = []

    # Use the generator method for coincidences
    # Choose which way to choose photon pairs within coincidences
    if PairMode == "Exclusive":
        for event in Generator:

            # Check if there are exactly two photons giving an LoR in acceptance
            if sd.TwoHitEvent( event, DetectorRadius, ZMin, ZMax ):
                hitRadius = sd.FindHitRadius( event, DetectorRadius )
                hitRadii.append( hitRadius )

                # Check truth information to distinguish randoms from scatters
                if sd.SameEventID( event ):
                    hitRadiiTruth.append( hitRadius )

    elif PairMode == "TakeAllGoods":
        for event in Generator:

            # Pair every subsequent photon with the first
            if len( event ) > 1:
                firstPhoton = event[0]
                for secondPhotonIndex in range( 1, len( event ) ):
                    #pair = np.take_along_axis( event, np.array([ 0, secondPhotonIndex ]), axis=0 )
                    pair = [ firstPhoton, event[secondPhotonIndex] ]
                    
                    # Now just repeat the "Exclusive" calculation
                    if sd.TwoHitEvent( pair, DetectorRadius, ZMin, ZMax ):
                        hitRadius = sd.FindHitRadius( pair, DetectorRadius )
                        hitRadii.append( hitRadius )
                        if sd.SameEventID( pair ):
                            hitRadiiTruth.append( hitRadius )
    else:
        print( "Unrecognised pairing mode: " + PairMode )
        return

    # First find the true and R+S coincidences
    y, x, patches = plt.hist( hitRadii, bins=26, range=[-130, 130] )
    plt.clf()
    NECRInSimWin, truesInSimWin, rPlusSInSimWin = NECRFromHistogram( x, y, SimulationWindow )
    
    # Then, find the number of scattered events by re-running NECRFromHistograms only for photons that have the same EventID
    y, x, patches = plt.hist( hitRadiiTruth, bins=26, range=[-130, 130] )
    plt.clf()
    _, _, scattersInSimWin = NECRFromHistogram( x, y, SimulationWindow )
    randomsInSimWin = rPlusSInSimWin - scattersInSimWin

    return NECRInSimWin, truesInSimWin, scattersInSimWin, randomsInSimWin
