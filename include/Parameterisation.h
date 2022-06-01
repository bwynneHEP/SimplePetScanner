#ifndef Parameterisation_h
#define Parameterisation_h 1

#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"

class Parameterisation : public G4VPVParameterisation
{
  public:
    Parameterisation( G4int DetectorsPerRing );
    ~Parameterisation(){};

    void ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const override;

  private:
    G4int m_detectorsPerRing;
};

#endif
