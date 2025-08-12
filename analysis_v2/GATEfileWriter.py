import SimulationDataset as sd
import numpy as np
import uproot
import math

# Data structure conversion
def PhotonToUpROOT( Photon, Suffix, ModuleIDs, BatchCounter, UpROOTdict ):

    UpROOTdict[ "time" + Suffix ][ BatchCounter ] = Photon[ sd.DATASET_TIME ]
    UpROOTdict[ "eventID" + Suffix ][ BatchCounter ] = Photon[ sd.DATASET_EVENT ]
    UpROOTdict[ "energy" + Suffix ][ BatchCounter ] = Photon[ sd.DATASET_ENERGY ]
    UpROOTdict[ "comptonPhantom" + Suffix ][ BatchCounter ] = 0 # TODO
    
    # These are helpful for debugging but not necessary
    UpROOTdict[ "globalPosX" + Suffix ][ BatchCounter ] = Photon[ sd.DATASET_R ] * math.cos( Photon[ sd.DATASET_PHI ] )
    UpROOTdict[ "globalPosY" + Suffix ][ BatchCounter ] = Photon[ sd.DATASET_R ] * math.sin( Photon[ sd.DATASET_PHI ] )
    UpROOTdict[ "globalPosZ" + Suffix ][ BatchCounter ] = Photon[ sd.DATASET_Z ]

    # Retrieve the module IDs
    # Note that it's fine to use a tuple () as a dict key, but not a list []
    globalCoordinates = ( Photon[ sd.DATASET_R ], Photon[ sd.DATASET_PHI ], Photon[ sd.DATASET_Z ] )
    if globalCoordinates not in ModuleIDs:
      #print( "Could not find module IDs for " + str( globalCoordinates ) )
      return

    # Only crystal-mode data contains full info
    photonModules = ModuleIDs[ globalCoordinates ]
    if len( photonModules ) < 3:
      #print( "Lacking full module ID info (" + str( len( photonModules ) ) + " fields found of 3 required)" )
      return

    # Save module IDs
    UpROOTdict[ "crystalID" + Suffix ][ BatchCounter ] = photonModules[ 2 ]   # "InBlock" in our terms
    UpROOTdict[ "submoduleID" + Suffix ][ BatchCounter ] = photonModules[ 0 ] # "Ring" in our terms
    UpROOTdict[ "moduleID" + Suffix ][ BatchCounter ] = 1 # Not needed?
    UpROOTdict[ "rsectorID" + Suffix ][ BatchCounter ] = photonModules[ 1 ]   # "Block" in our terms


# Run a coincidence generator and write results to a file
def GATEfromGenerator( BatchSize, Generator, FileName, ModuleIDs, DetectorRadius, PairMode="Exclusive", ZMin=0.0, ZMax=0.0 ):

    # Check for nonsense batch size
    if ( BatchSize < 16 ):
        print( "Requested BatchSize " + str( BatchSize ) + " is too small: using 256" )
        BatchSize = 256

    # Create output file
    outputFile = uproot.recreate( FileName )
    outputFile.mktree( "Coincidences", { "time1": np.float64,
                                         "time2": np.float64,
                                         "eventID1": np.int32,
                                         "eventID2": np.int32,
                                         "energy1": np.float32,
                                         "energy2": np.float32,
                                         "globalPosX1": np.float32, # debug
                                         "globalPosX2": np.float32, # debug
                                         "globalPosY1": np.float32, # debug
                                         "globalPosY2": np.float32, # debug
                                         "globalPosZ1": np.float32, # debug
                                         "globalPosZ2": np.float32, # debug
                                         "comptonPhantom1": np.int32,
                                         "comptonPhantom2": np.int32,
                                         "crystalID1": np.int32,
                                         "crystalID2": np.int32,
                                         "submoduleID1": np.int32,
                                         "submoduleID2": np.int32,
                                         "moduleID1": np.int32,
                                         "moduleID2": np.int32,
                                         "rsectorID1": np.int32,
                                         "rsectorID2": np.int32 })

    # Initialise output structure
    outputDict = { "time1": np.zeros( BatchSize, np.float64 ),
                   "time2": np.zeros( BatchSize, np.float64 ),
                   "eventID1": np.zeros( BatchSize, np.int32 ),
                   "eventID2": np.zeros( BatchSize, np.int32 ),
                   "energy1": np.zeros( BatchSize, np.float32 ),
                   "energy2": np.zeros( BatchSize, np.float32 ),
                   "globalPosX1": np.zeros( BatchSize, np.float32 ), # debug
                   "globalPosX2": np.zeros( BatchSize, np.float32 ), # debug
                   "globalPosY1": np.zeros( BatchSize, np.float32 ), # debug
                   "globalPosY2": np.zeros( BatchSize, np.float32 ), # debug
                   "globalPosZ1": np.zeros( BatchSize, np.float32 ), # debug
                   "globalPosZ2": np.zeros( BatchSize, np.float32 ), # debug
                   "comptonPhantom1": np.zeros( BatchSize, np.int32 ),
                   "comptonPhantom2": np.zeros( BatchSize, np.int32 ),
                   "crystalID1": np.zeros( BatchSize, np.int32 ),
                   "crystalID2": np.zeros( BatchSize, np.int32 ),
                   "submoduleID1": np.zeros( BatchSize, np.int32 ),
                   "submoduleID2": np.zeros( BatchSize, np.int32 ),
                   "moduleID1": np.zeros( BatchSize, np.int32 ),
                   "moduleID2": np.zeros( BatchSize, np.int32 ),
                   "rsectorID1": np.zeros( BatchSize, np.int32 ),
                   "rsectorID2": np.zeros( BatchSize, np.int32 ) }

    # Keeping track of place within batch
    batchCounter = 0

    # Use the generator method for coincidences
    # Choose which way to choose photon pairs within coincidences
    if PairMode == "Exclusive":
        for event in Generator:

            # Check if there are exactly two photons giving an LoR in acceptance
            if sd.TwoHitEvent( event, DetectorRadius, ZMin, ZMax ):
                PhotonToUpROOT( event[0], "1", ModuleIDs, batchCounter, outputDict )
                PhotonToUpROOT( event[1], "2", ModuleIDs, batchCounter, outputDict )
                batchCounter += 1

               # Check for end of batch
                if batchCounter == BatchSize:
                    outputFile[ "Coincidences" ].extend( outputDict )
                    batchCounter = 0

        # Tidy up the last entries
        for key in outputDict:
          outputDict[ key ] = outputDict[ key ][ 0:batchCounter ]
        outputFile[ "Coincidences" ].extend( outputDict )

    elif PairMode == "TakeAllGoods":
        for event in Generator:

            # Pair every subsequent photon with the first
            if len( event ) > 1:
                firstPhoton = event[0]
                for secondPhotonIndex in range( 1, len( event ) ):
                    pair = [ firstPhoton, event[secondPhotonIndex] ]
                    
                    # Now just repeat the "Exclusive" calculation
                    if sd.TwoHitEvent( pair, DetectorRadius, ZMin, ZMax ):
                        PhotonToUpROOT( pair[0], "1", ModuleIDs, batchCounter, outputDict )
                        PhotonToUpROOT( pair[1], "2", ModuleIDs, batchCounter, outputDict )
                        batchCounter += 1

                        # Check for end of batch
                        if batchCounter == BatchSize:
                            outputFile[ "Coincidences" ].extend( outputDict )
                            batchCounter = 0

        # Tidy up the last entries
        for key in outputDict:
          outputDict[ key ] = outputDict[ key ][ 0:batchCounter ]
        outputFile[ "Coincidences" ].extend( outputDict )

    else:
        print( "Unrecognised pairing mode: " + PairMode )
        return

