#include "ActionInitialization.h"
#include "GeneratorAction.h"
#include "RunAction.h"
#include "EventAction.h"

ActionInitialization::ActionInitialization() : G4VUserActionInitialization()
{
}

ActionInitialization::~ActionInitialization()
{
}

// Three actions to set up
// - generating particles
// - controlling the whole run
// - controlling a single event
void ActionInitialization::Build() const
{
  this->SetUserAction( new GeneratorAction() );
  this->SetUserAction( new RunAction() );
  this->SetUserAction( new EventAction() );
}
