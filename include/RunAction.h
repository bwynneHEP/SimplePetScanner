#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    ~RunAction() override;

    void BeginOfRunAction( const G4Run* ) override;
    void EndOfRunAction( const G4Run* ) override;
};

#endif
