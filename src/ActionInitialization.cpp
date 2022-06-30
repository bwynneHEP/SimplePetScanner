#include "ActionInitialization.h"
#include "LinearSourceAction.h"

ActionInitialization::ActionInitialization( DecayTimeFinderAction * decayTimeFinder ) : G4VUserActionInitialization(),
  m_decayTimeFinder( decayTimeFinder )
{
}

ActionInitialization::~ActionInitialization()
{
}

void ActionInitialization::Build() const
{
  this->SetUserAction( new LinearSourceAction( -350.0, 350.0 ) );
  this->SetUserAction( m_decayTimeFinder );
}
