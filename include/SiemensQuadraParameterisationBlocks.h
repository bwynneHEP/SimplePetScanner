// This class describes how the detector crystals are laid out in space

#ifndef SiemensQuadraParameterisationBlocks_h
#define SiemensQuadraParameterisationBlocks_h 1

#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"

class SiemensQuadraParameterisationBlocks : public G4VPVParameterisation
{
  public:
    SiemensQuadraParameterisationBlocks();
    ~SiemensQuadraParameterisationBlocks(){};

    void ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const override;
};

#endif
