#include "ActionInitialization.h"
#include "GeneratorAction.h"

ActionInitialization::ActionInitialization() : G4VUserActionInitialization()
{
}

ActionInitialization::~ActionInitialization()
{
}

void ActionInitialization::Build() const
{
  this->SetUserAction( new GeneratorAction() );
}
