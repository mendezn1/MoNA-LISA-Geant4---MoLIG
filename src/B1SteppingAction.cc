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
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>
#include <vector>

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) {
    const B1DetectorConstruction* detectorConstruction
      = static_cast<const B1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();
  }

  // get volume of the current step
  G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

  // check if we are in scoring volume
  //if (volume != fScoringVolume) return;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);

  G4ThreeVector ex_momentum;
  G4double ex_energy;

  G4ThreeVector fex_momentum;
  G4double fex_energy;
  G4ThreeVector nex_momentum;
  G4double nex_energy;
  G4double fex_mass;
  G4double nex_mass;
  G4String particle_name;
  vector<double> neutron_energy;
  vector<double> neutron_momentum_x; //momentum-x array
  vector<double> neutron_momentum_y; //momentum-y array
  vector<double> neutron_momentum_z; //momentum-z array
  vector<double> fragment_energy;
  vector<double> fragment_momentum_x; //momentum-x array
  vector<double> fragment_momentum_y; //momentum-y array
  vector<double> fragment_momentum_z; //momentum-z array

  //ex_momentum = step->GetPostStepPoint()->GetMomentum();
  //ex_energy   = step->GetPostStepPoint()->GetTotalEnergy();

  particle_name = step->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
  if (particle_name == "neutron")
  {
    //nex_mass = step->GetTrack()->GetParticleDefinition()->GetPDGMass();
    G4double energy_n_tot = step->GetPostStepPoint()->GetTotalEnergy();
    G4double n_momentum_x = step->GetPostStepPoint()->GetMomentum().x();
    G4double n_momentum_y = step->GetPostStepPoint()->GetMomentum().y();
    G4double n_momentum_z = step->GetPostStepPoint()->GetMomentum().z();
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
   G4double energy_f_tot = step->GetPostStepPoint()->GetTotalEnergy();
   G4double f_momentum_x = step->GetPostStepPoint()->GetMomentum().x();
   G4double f_momentum_y = step->GetPostStepPoint()->GetMomentum().y();
   G4double f_momentum_z = step->GetPostStepPoint()->GetMomentum().z();
   //fex_momentum = ex_momentum;
   //fex_energy = ex_energy;
   fragment_energy.push_back(energy_f_tot);
   fragment_momentum_x.push_back(f_momentum_x);
   fragment_momentum_y.push_back(f_momentum_y);
   fragment_momentum_z.push_back(f_momentum_z);
  }

  //G4cout << particle_name <<  G4endl;
  //G4cout << "LAKSJDFHALSKDJFHALSKDJFHALSKDJFHASLKDJFHASLDKFJHASLDKFJHASDLFKJSHDAFLKASJF" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
