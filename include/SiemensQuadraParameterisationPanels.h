// This class describes how the detector crystals are laid out in space

#ifndef SiemensQuadraParameterisationPanels_h
#define SiemensQuadraParameterisationPanels_h 1

#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"

class SiemensQuadraParameterisationPanels : public G4VPVParameterisation
{
  public:
    SiemensQuadraParameterisationPanels();
    ~SiemensQuadraParameterisationPanels(){};

    void ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const override;
};

#endif
