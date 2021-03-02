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

#include "BasicEventAction.hh"
#include "BasicRunAction.hh"
#include "BasicPETSD.hh"
#include "BasicPETHit.hh"
#include "BasicAnalysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "G4THitsMap.hh"

#include "Randomize.hh"
#include <iomanip>

//

BasicEventAction::BasicEventAction(BasicRunAction* runAction)
 : G4UserEventAction(),
   fRunAction(runAction),
   fDetHCID(-1),
   fPhanHCID(-1),
   fSumEdep(0.)
{ }

//

BasicEventAction::~BasicEventAction()
{}

//

BasicPETHitsCollection*
BasicEventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection
    = static_cast<BasicPETHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));

  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID;
    G4Exception("BasicEventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }

  return hitsCollection;
}

//

void BasicEventAction::PrintEventStatistics(
                              G4double detectorEdep, G4double detectorTrackLength) const
{
  // print event statistics
  G4cout
     << "   Detector: total energy: "
     << std::setw(7) << G4BestUnit(detectorEdep, "Energy")
     << "       total track length: "
     << std::setw(7) << G4BestUnit(detectorTrackLength, "Length")
     << G4endl;
}

//

void BasicEventAction::BeginOfEventAction(const G4Event* /*event*/)
{}

//

void BasicEventAction::EndOfEventAction(const G4Event* event)
{
  // Get hits collections IDs (only once)
  if ( fDetHCID == -1 ) {
    fDetHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("DetectorHitsCollection");
    G4cout << "\n fDetHCID: " << fDetHCID << "\n" << G4endl;   
  }

  // Get hits collections IDs (only once)
  if ( fPhanHCID < 0 ) {
    fPhanHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("patient/edep");
    G4cout << "\n fPhanHCID: " << fPhanHCID << "\n" << G4endl;   
  }
  
  G4int evtNb = event->GetEventID();
  
  G4cout << G4endl << "evtNb = " << evtNb ;

  // Get hits collections
  auto detHC = GetHitsCollection(fDetHCID, event);

  // Get hit with total values
  auto detHit = (*detHC)[detHC->entries()-1];

  // get deposited energy
  G4double dep = detHit->GetEdep();

  // I have followed the code in B3bRun.ccc
  //  
  G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  if(!HCE) return;
  
  G4THitsMap<G4double>* evtMap = 
    static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fPhanHCID));
  
  std::map<G4int,G4double*>::iterator itr;

  G4double phanedep = 0.;
  G4double meanedep = 0.;
  G4int itrs = 0;
  
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    phanedep = *(itr->second);
    
    G4int copyNb  = (itr->first);
    //if (edep > 10.)
    G4cout.precision(8);
    G4cout << G4endl << "  patient " << copyNb << ": " << phanedep/keV << " keV " << G4endl << " sum edep: " << fSumEdep;
        
        // get analysis manager 
        auto analysisManager = G4AnalysisManager::Instance();

        // fill phantom histogram data
        analysisManager->FillH1(0, phanedep);
        analysisManager->FillNtupleDColumn(0, phanedep);
        analysisManager->AddNtupleRow();

    fSumEdep += phanedep/keV; 
    itrs++; 
  }  

  meanedep = fSumEdep/(evtNb + 1);

  G4cout.precision(8); 
  G4cout << " Mean Edep: " << meanedep << G4endl;
  
  // defining a Good Event
  G4double EnergyRes = 1.022*0.2;
  G4double Threshold = (1.022 - EnergyRes)*MeV;
  G4double ScatterThreshold = 10.0*keV;
  if (dep > Threshold) fRunAction->CountEvent();  
  if (phanedep > ScatterThreshold) fRunAction->ScatterEvent();

  // Print per event (modulo n)
  auto eventID = event->GetEventID();
  auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
  if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {
    G4cout << "---> End of event: " << eventID << G4endl;

    PrintEventStatistics(
      dep, detHit->GetTrackLength());
  }

  // get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // fill histograms
  analysisManager->FillH1(0, dep);
  analysisManager->FillH1(1, detHit->GetTrackLength());
  //analysisManager->FillH1(0, phanedep);

  analysisManager->FillNtupleDColumn(0, dep);
  analysisManager->FillNtupleDColumn(1, detHit->GetTrackLength());
  analysisManager->AddNtupleRow();

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...


void BasicEventAction::Merge(const G4UserEventAction* eventAction)
{
  const BasicEventAction* localRun = static_cast<const BasicEventAction*>(eventAction);
  fSumEdep += localRun->fSumEdep;

  BasicEventAction::Merge(eventAction);
} 





// energy resolution for LSO is 10.6%
// good event = edep > 89.4% of 1.022MeV
