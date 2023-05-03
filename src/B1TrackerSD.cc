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
//
/// \file B1TrackerSD.cc
/// \brief Implementation of the B1TrackerSD class

#include "B1TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"


#include <vector> // included by NM 2023-3-20
#include "edecay_calc_SD.hh" // included by NM 2023-3-20
#include "HistoManager.hh"

using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1TrackerSD::B1TrackerSD(const G4String& name,
                         const G4String& hitsCollectionName)
 : G4VSensitiveDetector(name),
   fHitsCollection(NULL)
{
  collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1TrackerSD::~B1TrackerSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1TrackerSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection

  fHitsCollection
    = new B1TrackerHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce

  G4int hcID
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double energy_n_tot;
G4double n_momentum_x;
G4double n_momentum_y;
G4double n_momentum_z;
G4double energy_f_tot;
G4double f_momentum_x;
G4double f_momentum_y;
G4double f_momentum_z;
vector<double> energy_tot;
vector<double> neutron_energy;
vector<double> neutron_momentum_x; //momentum-x array
vector<double> neutron_momentum_y; //momentum-y array
vector<double> neutron_momentum_z; //momentum-z array
vector<double> fragment_energy;
vector<double> fragment_momentum_x; //momentum-x array
vector<double> fragment_momentum_y; //momentum-y array
vector<double> fragment_momentum_z; //momentum-z array

G4bool B1TrackerSD::ProcessHits(G4Step* aStep,
                                     G4TouchableHistory*)
{
  // Predefined quantities
  G4double edep = aStep->GetPostStepPoint()->GetTotalEnergy();
  G4String particle_name = aStep->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();

  if (edep==0.) return false;

  B1TrackerHit* newHit = new B1TrackerHit();

  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  newHit->SetTrackName(particle_name);
  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()
                                               ->GetCopyNumber());
  newHit->SetEdep(edep);
  newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());
  newHit->SetMom (aStep->GetPostStepPoint()->GetMomentum());

  fHitsCollection->insert( newHit );

  //newHit->Print();

  if (particle_name == "neutron")
  {
    //nex_mass = step->GetTrack()->GetParticleDefinition()->GetPDGMass();
    energy_n_tot = aStep->GetPostStepPoint()->GetTotalEnergy();
    n_momentum_x = aStep->GetPostStepPoint()->GetMomentum().x();
    n_momentum_y = aStep->GetPostStepPoint()->GetMomentum().y();
    n_momentum_z = aStep->GetPostStepPoint()->GetMomentum().z();
    //nex_energy = ex_energy;
    //fEventAction->StoreNeutronEnergy(energy_n_tot);
    neutron_energy.push_back(energy_n_tot);
    neutron_momentum_x.push_back(n_momentum_x);
    neutron_momentum_y.push_back(n_momentum_y);
    neutron_momentum_z.push_back(n_momentum_z);
  }

  if (particle_name == "Ca52")
  {
   //fex_mass = step->GetTrack()->GetParticleDefinition()->GetPDGMass();
   energy_f_tot = aStep->GetPostStepPoint()->GetTotalEnergy();
   f_momentum_x = aStep->GetPostStepPoint()->GetMomentum().x();
   f_momentum_y = aStep->GetPostStepPoint()->GetMomentum().y();
   f_momentum_z = aStep->GetPostStepPoint()->GetMomentum().z();
   //fex_momentum = ex_momentum;
   //fex_energy = ex_energy;
   fragment_energy.push_back(energy_f_tot);
   fragment_momentum_x.push_back(f_momentum_x);
   fragment_momentum_y.push_back(f_momentum_y);
   fragment_momentum_z.push_back(f_momentum_z);
  }

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void B1TrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( true ) {
     G4int nofHits = fHitsCollection->entries();
     G4cout << G4endl
            << "-------->Hits Collection: in this event they are " << nofHits
            << " hits in the tracker chambers: " << G4endl;
     for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
     //for ( G4int i=0; i<vnts; i++ ) {
     // G4cout << neutron_energy[i]*fragment_energy[i] << G4endl;
     //}
     int vnts = 9e3;
     //test
     vector<double> edecays_geant = edecay_v2(vnts,neutron_energy ,neutron_momentum_x, neutron_momentum_y, neutron_momentum_z,
                                              fragment_energy,fragment_momentum_x,fragment_momentum_y,fragment_momentum_z);
     for (auto& edecay_vals : edecays_geant){
      G4AnalysisManager::Instance()->FillH1(1,edecay_vals);
     }
     for (auto& neutron_g4_Etot_vals : neutron_energy){
      G4AnalysisManager::Instance()->FillH1(3,neutron_g4_Etot_vals);
     }
     for (auto& fragment_g4_Etot_vals : fragment_energy){
      G4AnalysisManager::Instance()->FillH1(5,fragment_g4_Etot_vals);
     }
     for (auto& neutron_g4_px_vals : neutron_momentum_x){
      G4AnalysisManager::Instance()->FillH1(7,neutron_g4_px_vals);
     }
     for (auto& neutron_g4_py_vals : neutron_momentum_y){
      G4AnalysisManager::Instance()->FillH1(9,neutron_g4_py_vals);
     }
     for (auto& neutron_g4_pz_vals : neutron_momentum_z){
      G4AnalysisManager::Instance()->FillH1(11,neutron_g4_pz_vals);
     }
     for (auto& fragment_g4_px_vals : fragment_momentum_x){
      G4AnalysisManager::Instance()->FillH1(13,fragment_g4_px_vals);
     }
     for (auto& fragment_g4_py_vals : fragment_momentum_y){
      G4AnalysisManager::Instance()->FillH1(15,fragment_g4_py_vals);
     }
     for (auto& fragment_g4_pz_vals : fragment_momentum_z){
      G4AnalysisManager::Instance()->FillH1(17,fragment_g4_pz_vals);
     }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
