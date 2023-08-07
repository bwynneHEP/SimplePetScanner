import sys
sys.path.append('../analysis/')
from SimulationDataset import *

import matplotlib.pyplot as mpl
params = {'legend.fontsize': 15,
          'legend.title_fontsize': 15,
          'legend.loc': "lower right",
          'axes.labelsize': 15,
          'xtick.labelsize': 15,
          'ytick.labelsize': 15}
mpl.rcParams.update(params)

def ReadBackground( crystalData ):
    
    z = []
    phi = []
    energyNearSide = []
    energyFarSide = []
    for dataset in crystalData:
        for i in range( dataset.size() ):
            event = dataset.SampleOneEvent()
            for hit in event:
                phiVal = hit[5]
            
                # Ensure it's not a decay in the same crystal
                if math.fabs( phiVal ) > math.pi / 4.0:
                    z.append( hit[6] )
                    energyFarSide.append( hit[2] )
                
                    # Rotate to opposite side centred
                    phiVal -= math.pi
                    if phiVal < -math.pi:
                        phiVal += 2.0 * math.pi
                    phi.append( phiVal )
                    
                else:
                    energyNearSide.append( hit[2] )
                    
    return energyNearSide, energyFarSide, z, phi

def Make2DHistForAttenuation( energy, z, phi ):

    print( "Detector-crossing hits: ", len(z) )
    hist = mpl.hist2d( z, phi, bins=[11,22], range=[[-550,550], [-math.pi,math.pi]], weights=energy )
    mpl.gcf().set_size_inches(10, 10)
    mpl.xlabel( "Axial coordinate [mm]" )
    mpl.ylabel( "Phi coordinate [radians]")
    mpl.show()

# Longish detector with a short PET phantom in the middle
datasetSize = 1000000
detectorLength = 1100
phantomLength = 100

detectorMaterial = "LSO"
crystalData = CreateDataset( detectorLength, "Siemens", phantomLength, "Siemens", datasetSize, 100.0, 900.0, detectorMaterial )
eNear, eFar, z, phi = ReadBackground( [crystalData] )
Make2DHistForAttenuation( eFar, z, phi )