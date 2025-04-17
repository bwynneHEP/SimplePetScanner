#include "DetectorConstruction.h"
#include "ActionInitialization.h"
#include "DecayTimeFinderAction.h"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4PhysListFactory.hh"
#include "LBE.hh"
// #include "QBBC.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "Randomize.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4HadronicParameters.hh"

void printHelp()
{
  std::cout << "Usage: build/SimplePetScanner [OPTIONS]\n"
              << "Options:\n"
              << "  --help                Show help message\n"
              << "  -n                    Number of events to simulate\n"
              << "  --gui                 Activate the gui\n"
              << "  --detector            Select detector geometry and granularity\n"
              << "  --detectorLengthMM    Set the detector length in mm\n"
              << "  --detectorMaterial    Set the scintillator material to use\n"
              << "  --source              Set the radioactive source (tracer or intrinsic detector background)\n"
              << "  --sourceOffsetMM      Shift the linear source by n mm\n"
              << "  --phantomLengthMM     Set the phantom length in mm\n"
              << "  --outputFileName      Override the default output file name\n"
              << "  --decayOutputFileName Create output containing the radioactive decay info\n"
              << "  --randomSeed          Override the default random seed\n"
              << "  --nAluminiumSleeves   Create sensitivity phantom with n sleeves\n"
              << "See README.md for more complete documentation\n";
} 

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
  G4double sourceOffsetMM = 0.0;
  G4int nAluminiumSleeves = 0;

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
    else if (argument == "--help") 
    { 
      printHelp();
      exit(0);
    }
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
    else if ( argument == "--sourceOffsetMM" )
    {
      if ( nextArgument.size() ) sourceOffsetMM = nextInteger;
      else
      {
        std::cerr << "Did not find source offset to use" << std::endl;
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
    else if ( argument == "--nAluminiumSleeves" )
    {
      if (nextArgument.size() ) {
        if (nextInteger > 0 && nextInteger <= 5)
          nAluminiumSleeves = nextInteger;
        else
        {
          std::cerr << "Invalid sleeve count, enter a value between 1 and 5" << std::endl;
          return 1;
        }
      }
      else 
      {
        std::cerr << "Did not find the number of aluminium sleeves to use" << std::endl;
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
  // G4VModularPhysicsList* physicsList = new QBBC();
  // Use QGSP_BIC_HP_EMY which contains G4EmStandardPhysics_option3 
  G4PhysListFactory physListFactory;
  G4VModularPhysicsList* physicsList = physListFactory.GetReferencePhysList("QGSP_BIC_HP_EMY");
  physicsList->RegisterPhysics( new G4StepLimiterPhysics() );
  // physicsList->RegisterPhysics( new G4RadioactiveDecayPhysics() ); // Already included in QGSP_BIC_HP_EMY, add only if using another physics list
  runManager->SetUserInitialization( physicsList );

  // Set user action classes
  DecayTimeFinderAction * decayTimeFinder = new DecayTimeFinderAction( decayOutputFileName );
  ActionInitialization * actions = new ActionInitialization( decayTimeFinder, sourceName, detectorLength, phantomLength, detectorMaterial, sourceOffsetMM );
  runManager->SetUserInitialization( actions );

  // Set up detector
  DetectorConstruction * detector = new DetectorConstruction( decayTimeFinder, detectorName, detectorLength, phantomLength, outputFileName, decayOutputFileName, detectorMaterial, nAluminiumSleeves, sourceOffsetMM );
  runManager->SetUserInitialization( detector );

  // Set up the macros
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  G4HadronicParameters::Instance()->SetTimeThresholdForRadioactiveDecay( 1.0e+60*CLHEP::year ); //increase the time threshold for ion decays
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
