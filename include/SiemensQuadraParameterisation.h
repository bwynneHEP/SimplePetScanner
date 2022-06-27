// This class describes how the detector crystals are laid out in space

#ifndef SiemensQuadraParameterisation_h
#define SiemensQuadraParameterisation_h 1

#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"

class SiemensQuadraParameterisation : public G4VPVParameterisation
{
  public:
    SiemensQuadraParameterisation();
    ~SiemensQuadraParameterisation(){};

    void ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const override;
};

#endif
