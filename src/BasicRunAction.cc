//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//

#include "BasicRunAction.hh"
#include "BasicAnalysis.hh"
#include "BasicDetectorConstruction.hh"
#include "BasicEventAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//


BasicRunAction::BasicRunAction()
 : G4UserRunAction()

{
  // set printing run number only
  G4RunManager::GetRunManager()->SetPrintProgress(0);

  // Create analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in BasicAnalysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories
  analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
  // Note: merging ntuples is available only with Root output

  // Creating histograms
  analysisManager->CreateH1("Energy","Energy Deposited", 50, 0.,1.25*MeV);
  analysisManager->CreateH1("Length","Track Length in Detector", 50, 0., 1.0*mm);
  analysisManager->CreateH1("Energy2","Energy Deposited in Phantom", 50, 0.,1.25*MeV);
  
  // Creating ntuple
  analysisManager->CreateNtuple("Basic", "Edep spacial distribution");
  analysisManager->CreateNtupleDColumn("Edep");
  analysisManager->CreateNtupleDColumn("TrackLength");
  analysisManager->CreateNtupleDColumn("PhanEdep");
  analysisManager->FinishNtuple();
}

//

BasicRunAction::~BasicRunAction()
{
  delete G4AnalysisManager::Instance();
}

//

void BasicRunAction::BeginOfRunAction(const G4Run* run)
{

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Reset the GoodEvent counter
  Reset();

  // Reset the Scatter Event counter
  ScatterReset();

  // Open an output file featuring the runID
  G4int runid = run->GetRunID();
  G4String fileName = "BasicOut" + std::to_string(runid);
  analysisManager->OpenFile(fileName);

}

//

void BasicRunAction::EndOfRunAction(const G4Run* run)
{

  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;



  const BasicEventAction* event = static_cast<const BasicEventAction*>(G4RunManager::GetRunManager()->GetUserEventAction());
  G4double sumEdep   = event->GetSumEdep(); 


  // print histogram statistics

  auto analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->GetH1(1) ) {

    G4int goodEvents = GoodEventCount;
    G4int scatEvents = ScatterEventCount;
    G4double sensitivity = (G4double(goodEvents)/nofEvents) * 100;
    G4double totaledep = (G4double(sumEdep/nofEvents)) ;
    G4double scatterfrac = (G4double(scatEvents)/nofEvents) * 100;

    G4cout << " Detector length: " << DetLength << " m " << G4endl;
    G4cout << " Crystal length: " << CrystLength << " cm " << G4endl;
    G4cout << " Good events: " << goodEvents << G4endl;
    G4cout << " Scatter events: " << scatEvents << G4endl;
    //G4cout << " Phantom Material: " << phantom_mat << G4endl;
    G4cout << " Crude sensitivity: " << std::setprecision(5) << sensitivity << " per cent " << G4endl;
    G4cout << " Crude Scatter Fraction: " << std::setprecision(5) << scatterfrac << " per cent " << G4endl;
    G4cout << " Mean Energy Deposited in Patient: " << totaledep << "keV " << G4endl;
  

    std::ofstream myfile;
    myfile.open("Data.csv", std::ofstream::app);
    //myfile << std::to_string(CrystLength)+", "+std::to_string(sensitivity) +"\n";
    myfile << std::to_string(DetLength)+", "+std::to_string(CrystLength)+", "+std::to_string(sensitivity)+", "+std::to_string(totaledep) +"\n";
    myfile.close();

  }


  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//
