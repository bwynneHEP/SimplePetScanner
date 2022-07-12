#include "ActionInitialization.h"
#include "LinearSourceAction.h"
#include "CrystalIntrinsicAction.h"

#include "G4SystemOfUnits.hh"

ActionInitialization::ActionInitialization( DecayTimeFinderAction * decayTimeFinder ) : G4VUserActionInitialization(),
  m_decayTimeFinder( decayTimeFinder )
{
}

ActionInitialization::~ActionInitialization()
{
}

void ActionInitialization::Build() const
{
  //this->SetUserAction( new LinearSourceAction( -350.0*mm, 350.0*mm ) );
  this->SetUserAction( new CrystalIntrinsicAction( -500.0*mm, 500.0*mm, 400.0*mm, 420.0*mm ) ); // 41cm radius is to the middle of the 2cm thick crystal - wrong?
  this->SetUserAction( m_decayTimeFinder );
}
