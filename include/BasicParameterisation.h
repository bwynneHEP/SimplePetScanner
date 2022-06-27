// This class describes how the detector crystals are laid out in space

#ifndef BasicParameterisation_h
#define BasicParameterisation_h 1

#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"

class BasicParameterisation : public G4VPVParameterisation
{
  public:
    BasicParameterisation( G4int DetectorsPerRing );
    ~BasicParameterisation(){};

    void ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const override;

  private:
    G4int m_detectorsPerRing;
};

#endif
