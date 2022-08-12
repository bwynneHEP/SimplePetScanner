// This class describes how the detector crystals are laid out in space

#ifndef ExplorerParameterisationPanels_h
#define ExplorerParameterisationPanels_h 1

#include "G4VPVParameterisation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"

class ExplorerParameterisationPanels : public G4VPVParameterisation
{
  public:
    ExplorerParameterisationPanels( G4int nCopies );
    ~ExplorerParameterisationPanels(){};

    void ComputeTransformation( const G4int copyNo, G4VPhysicalVolume* physVol ) const override;

  private:
    std::vector< G4ThreeVector > m_positions;
    std::vector< G4RotationMatrix* > m_rotations;
    std::vector< G4VisAttributes* > m_visions;
};

#endif
