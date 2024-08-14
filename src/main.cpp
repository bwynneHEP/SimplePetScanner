#include "DetectorConstruction.h"
#include "ActionInitialization.h"
#include "DecayTimeFinderAction.h"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "QGSP_BERT_HP.hh"
#include "LBE.hh"
#include "QBBC.hh"
#include "G4OpticalPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "Randomize.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"


int main( int argc, char* argv[] )
{
  // Parse command-line arguments
  G4bool useGUI = false;
  G4int nEvents = 0;
  G4int randomSeed = 1234;
  std::string detectorMaterial = "";
  std::string detectorName = "";
  std::string sourceName = "";
  std::string outputFileName = "hits.csv";
  std::string decayOutputFileName = "";
  G4double detectorLength = -1.0;
  G4double phantomLength = -1.0;
  for ( int argi = 1; argi < argc; ++argi )
  {
    std::string argument = argv[ argi ];

    // Find argument parameter
    std::string nextArgument = "";
    if ( argi + 1 < argc )
    {
      nextArgument = argv[ argi + 1 ];
      if ( nextArgument[0] == '-' ) nextArgument = "";
      else ++argi;
    }

    // Find integer argument
    int nextInteger = 0;
    if ( nextArgument.size() )
    {
      nextInteger = atoi( nextArgument.c_str() );
    }

    // Examine arguments
    if ( argument == "--gui" ) useGUI = true;
    else if ( argument == "-n" )
    {
      if ( nextInteger ) nEvents = nextInteger;
      else
      {
        std::cerr << "Failed to parse number of events" << std::endl;
        return 1;
      }
    }
    else if ( argument == "--detector" )
    {
      if ( nextArgument.size() ) detectorName = nextArgument;
      else
      {
        std::cerr << "Did not find detector name to use" << std::endl;
        return 1;
      }
    }
    else if ( argument == "--detectorLengthMM" )
    {
      if ( nextArgument.size() ) detectorLength = nextInteger;
      else
      {
        std::cerr << "Did not find detector length to use" << std::endl;
        return 1;
      }
    }
    else if ( argument == "--detectorMaterial" )
    {
      if ( nextArgument.size() ) detectorMaterial = nextArgument;
      else
      {
        std::cerr << "Did not find detector material to use" << std::endl;
        return 1;
      }
    }
    else if ( argument == "--source" )
    {
      if ( nextArgument.size() ) sourceName = nextArgument;
      else
      {
        std::cerr << "Did not find source name to use" << std::endl;
        return 1;
      }
    }
    else if ( argument == "--phantomLengthMM" )
    {
      if ( nextArgument.size() ) phantomLength = nextInteger;
      else
      {
        std::cerr << "Did not find phantom length to use" << std::endl;
        return 1;
      }
    }
    else if ( argument == "--outputFileName" )
    {
      if ( nextArgument.size() ) outputFileName = nextArgument;
      else
      {
        std::cerr << "Did not find output file name to use" << std::endl;
        return 1;
      }
    }
    else if ( argument == "--decayOutputFileName" )
    {
      if ( nextArgument.size() ) decayOutputFileName = nextArgument;
      else
      {
        std::cerr << "Did not find decay output file name to use" << std::endl;
        return 1;
      }
    }
    else if ( argument == "--randomSeed" )
    {
      if ( nextArgument.size() ) randomSeed = nextInteger;
      else
      {
        std::cerr << "Did not find random seed to use" << std::endl;
        return 1;
      }
    }
    else
    {
      std::cerr << "Unrecognised argument: " << argument << std::endl;
      return 1;
    }
  }

  // Set up the random number generator
  G4Random::setTheEngine( new CLHEP::RanecuEngine );
  G4Random::setTheSeed( randomSeed );

  // Construct the default run manager
  G4RunManager* runManager = new G4RunManager();

  // Set up physics processes
  // G4VModularPhysicsList* physicsList = new QGSP_BERT_HP();
  // G4VModularPhysicsList* physicsList = new LBE();
  // Hanna: Based on initial studies chose to use QBBC 
  G4VModularPhysicsList* physicsList = new QBBC();
  physicsList->RegisterPhysics( new G4StepLimiterPhysics() );
  physicsList->RegisterPhysics( new G4RadioactiveDecayPhysics() ); // For the tracers
  physicsList->RegisterPhysics( new G4OpticalPhysics() ); // For optical photons
  runManager->SetUserInitialization( physicsList );

  // Set user action classes
  DecayTimeFinderAction * decayTimeFinder = new DecayTimeFinderAction( decayOutputFileName );
  ActionInitialization * actions = new ActionInitialization( decayTimeFinder, sourceName, detectorLength, phantomLength, detectorMaterial );
  runManager->SetUserInitialization( actions );

  // Set up detector
  DetectorConstruction * detector = new DetectorConstruction( decayTimeFinder, detectorName, detectorLength, phantomLength, outputFileName, decayOutputFileName, detectorMaterial );
  runManager->SetUserInitialization( detector );

  // Set up the macros
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand( "/run/initialize" );
  UImanager->ApplyCommand( "/control/execute run.mac" );
  UImanager->ApplyCommand( "/run/beamOn " + std::to_string( nEvents ) ); // even if it's zero, useful to initialise physics

  // Set up GUI
  G4VisManager* visManager = nullptr;
  if ( useGUI )
  {
    G4UIExecutive* ui = new G4UIExecutive( argc, argv );
    visManager = new G4VisExecutive();
    visManager->Initialize();
    UImanager->ApplyCommand( "/control/execute vis.mac" );
    ui->SessionStart();
    delete ui;
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted
  // in the main() program !
  //
  if ( useGUI ) delete visManager;
  delete runManager;

  return 0;
}
