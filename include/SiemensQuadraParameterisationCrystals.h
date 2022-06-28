// This class describes how the detector crystals are laid out in space

#ifndef SiemensQuadraParameterisationCrystals_h
#define SiemensQuadraParameterisationCrystals_h 1

#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"

class SiemensQuadraParameterisationCrystals : public G4VPVParameterisation
{
  public:
    SiemensQuadraParameterisationCrystals();
    ~SiemensQuadraParameterisationCrystals(){};

    void ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const override;
};

#endif
