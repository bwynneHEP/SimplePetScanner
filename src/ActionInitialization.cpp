#include "ActionInitialization.h"
#include "LinearSourceAction.h"
#include "CrystalIntrinsicAction.h"
#include "ExplorerDetector.h"
#include "SiemensQuadraDetector.h"

#include "G4SystemOfUnits.hh"

ActionInitialization::ActionInitialization( DecayTimeFinderAction * decayTimeFinder, std::string sourceName, G4double detectorLength, G4double phantomLength, std::string detectorMaterial)
  : G4VUserActionInitialization()
  , m_decayTimeFinder( decayTimeFinder )
  , m_sourceName( sourceName )
  , m_detectorLength( detectorLength )
  , m_phantomLength( phantomLength )
  , m_detectorMaterial( detectorMaterial )
{
}

ActionInitialization::~ActionInitialization()
{
}

void ActionInitialization::Build() const
{
  // std::cout << "Source name: " << m_sourceName << std::endl;
  if ( m_sourceName.substr( 0, 6 ) == "Linear" )
  {
    G4double phantomLength = 350.0*mm;
    if ( m_phantomLength > 0.0 )
    {
      phantomLength = m_phantomLength * mm / 2.0; // half-lengths
    }
    else if ( m_phantomLength == 0.0 )
    {
      std::cerr << "Cannot use a zero-length phantom as a source" << std::endl;
      exit(1);
    }
    this->SetUserAction( new LinearSourceAction( -phantomLength, phantomLength, m_sourceName.substr( 6 ) ) );
  }
  else if ( m_sourceName == "Siemens" )
  {
    G4double detectorLength = 512.0*mm;
    if ( m_detectorLength > 0.0 )
    {
      detectorLength = SiemensQuadraDetector::LengthForNRings( SiemensQuadraDetector::NRingsInLength( m_detectorLength ) ); // discrete length steps given by rings
      detectorLength /= 2.0; // half-lengths
    }
    if (m_detectorMaterial == "NaI")
      this->SetUserAction( new CrystalIntrinsicAction( -detectorLength, detectorLength, "NaI", 400.0*mm, 420.0*mm ) );
    else if (m_detectorMaterial == "BGO")
      this->SetUserAction( new CrystalIntrinsicAction( -detectorLength, detectorLength, "BGO", 400.0*mm, 420.0*mm ) );
    else if (m_detectorMaterial == "CsF")
      this->SetUserAction( new CrystalIntrinsicAction( -detectorLength, detectorLength, "CsF", 400.0*mm, 420.0*mm ) );
    else 
      this->SetUserAction( new CrystalIntrinsicAction( -detectorLength, detectorLength, "LSO", 400.0*mm, 420.0*mm ) );  
  }
  else if ( m_sourceName == "Explorer" )
  {
    G4double detectorLength = 936.46*mm;
    if ( m_detectorLength > 0.0 )
    {
      detectorLength = ExplorerDetector::LengthForNRings( ExplorerDetector::NRingsInLength( m_detectorLength ) ); // discrete length steps given by rings
      detectorLength /= 2.0; // half-lengths
    }
    if (m_detectorMaterial == "NaI")
      this->SetUserAction( new CrystalIntrinsicAction( -detectorLength, detectorLength, "NaI", 393*mm, 411.1*mm ) );
    else if (m_detectorMaterial == "BGO")
      this->SetUserAction( new CrystalIntrinsicAction( -detectorLength, detectorLength, "BGO", 393*mm, 411.1*mm ) );
    else if (m_detectorMaterial == "CsF")
      this->SetUserAction( new CrystalIntrinsicAction( -detectorLength, detectorLength, "CsF", 393*mm, 411.1*mm ) );
    else 
      this->SetUserAction( new CrystalIntrinsicAction( -detectorLength, detectorLength, "LYSO", 393*mm, 411.1*mm ) );
  }
  else
  {
    std::cerr << "Unrecognised source name: " << m_sourceName << std::endl;
    exit(1);
  }

  this->SetUserAction( m_decayTimeFinder );
}
