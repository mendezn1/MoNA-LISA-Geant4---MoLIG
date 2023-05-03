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
/// \file B1PrimaryGeneratorAction.cc
/// \brief Implementation of the B1PrimaryGeneratorAction class

//#include <TH1D.h>

#include "B1PrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "B1SteppingAction.hh"

#include "G4IonTable.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include "kinematics.hh"
#include "HistoManager.hh"
#include "edecay_calc.hh"
#include "B1EventAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

B1PrimaryGeneratorAction::B1PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fEnvelopeBox(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1PrimaryGeneratorAction::~B1PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  //this function is called at the begining of ecah event
  //

  // In order to avoid dependence of PrimaryGeneratorAction
  // on DetectorConstruction class we get Envelope volume
  // from G4LogicalVolumeStore.
  int events = 3e3;

  vector<vector<vector<double>>> data_92 = kinematics(events, 0.5, 0.5, 4, 52);
  vector<vector<vector<double>>> data_52 = kinematics(events, 1.9, 1.8, 2, 52);
  vector<vector<vector<double>>> data_32 = kinematics(events, 2.5, 3.6, 2, 52);
  vector<vector<double>> neutrons_92  = data_92[0];
  vector<vector<double>> fragments_92 = data_92[1];
  vector<vector<double>> neutrons_52  = data_52[0];
  vector<vector<double>> fragments_52 = data_52[1];
  vector<vector<double>> neutrons_32  = data_32[0];
  vector<vector<double>> fragments_32 = data_32[1];

  //string line;                    /* string to hold each line */
  //vector<vector<double>> array;    /* vector of vector<float> for 2d array */

  //string fragment_data = "boosted_data_fragments.txt";

  //ifstream f (fragment_data);   /* open file */
  //  if (!f.is_open()) {     /* validate file open for reading */
  //      perror (("error while opening file " + string(fragment_data)).c_str());
  //      return;
  //}

  //while (getline (f, line)) {         /* read each line */
  //    string val;                     /* string to hold value */
  //    vector<double> row;              /* vector for row of values */
  //    stringstream ss (line);          /* stringstream to parse csv */
  //    while (getline (ss, val, ','))   /* for each value */
  //        row.push_back (stod(val));  /* convert to float, add to row */
  //    array.push_back (row);          /* add row to array */
  //}
  //f.close();

  vector<double> neutron_energy_92; //energy array
  vector<double> neutron_px_92; //momentum-x array
  vector<double> neutron_py_92; //momentum-y array
  vector<double> neutron_pz_92; //momentum-z array
  vector<double> neutron_energy_52; //energy array
  vector<double> neutron_px_52; //momentum-x array
  vector<double> neutron_py_52; //momentum-y array
  vector<double> neutron_pz_52; //momentum-z array
  vector<double> neutron_energy_32; //energy array
  vector<double> neutron_px_32; //momentum-x array
  vector<double> neutron_py_32; //momentum-y array
  vector<double> neutron_pz_32; //momentum-z array

  vector<double> fragment_energy_92; //energy array
  vector<double> fragment_px_92; //momentum-x array
  vector<double> fragment_py_92; //momentum-y array
  vector<double> fragment_pz_92; //momentum-z array
  vector<double> fragment_energy_52; //energy array
  vector<double> fragment_px_52; //momentum-x array
  vector<double> fragment_py_52; //momentum-y array
  vector<double> fragment_pz_52; //momentum-z array
  vector<double> fragment_energy_32; //energy array
  vector<double> fragment_px_32; //momentum-x array
  vector<double> fragment_py_32; //momentum-y array
  vector<double> fragment_pz_32; //momentum-z array

  for (auto& row : fragments_92) {         /* read each line */
    fragment_energy_92.push_back(row[3]);          /* add row to array */
    fragment_px_92.push_back(row[0]);
    fragment_py_92.push_back(row[1]);
    fragment_pz_92.push_back(row[2]);
  }
  for (auto& row : fragments_52) {         /* read each line */
    fragment_energy_52.push_back(row[3]);          /* add row to array */
    fragment_px_52.push_back(row[0]);
    fragment_py_52.push_back(row[1]);
    fragment_pz_52.push_back(row[2]);
  }
  for (auto& row : fragments_32) {         /* read each line */
    fragment_energy_32.push_back(row[3]);          /* add row to array */
    fragment_px_32.push_back(row[0]);
    fragment_py_32.push_back(row[1]);
    fragment_pz_32.push_back(row[2]);
  }
  for (auto& row : neutrons_92) {         /* read each line */
    neutron_energy_92.push_back(row[3]);          /* add row to array */
    neutron_px_92.push_back(row[0]);
    neutron_py_92.push_back(row[1]);
    neutron_pz_92.push_back(row[2]);
  }
  for (auto& row : neutrons_52) {         /* read each line */
    neutron_energy_52.push_back(row[3]);          /* add row to array */
    neutron_px_52.push_back(row[0]);
    neutron_py_52.push_back(row[1]);
    neutron_pz_52.push_back(row[2]);
  }
  for (auto& row : neutrons_32) {         /* read each line */
    neutron_energy_32.push_back(row[3]);          /* add row to array */
    neutron_px_32.push_back(row[0]);
    neutron_py_32.push_back(row[1]);
    neutron_pz_32.push_back(row[2]);
  }

  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";

  cout << "\n";
  cout << "\n";
  cout << "\n";
  cout << "\n";
  cout << "\n";

  cout << "OH MY GOSH I LOVE GEANT IT IS THE MOST USEFUL AND EASIEST TOOL TO USE I DONT KNOW WHY PEOPLE DONT LIKE THIS TOOL" << G4endl;
  //for (auto& vals : neutron_energy_92){
  //  G4cout << vals << G4endl;
  //}
  //for (auto& vals : neutron_px_92){
  //  G4cout << vals << G4endl;
  //}
  //for (auto& vals : neutron_py_92){
  //  G4cout << vals << G4endl;
  //}
  //for (auto& vals : neutron_energy_92){
  //  G4cout << fixed << setprecision(3) << vals << G4endl;
  //}

  cout << "\n";
  cout << "\n";
  cout << "\n";
  cout << "\n";
  cout << "\n";

  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";
  cout << "da;sdkfja;sdfkljasd;flkajsdf;alksdjfa;ldkjfa;sdklfja;sdklfjas;dlkfja;sdlkfja;sdklfjas;dklfj";

  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  //Defining Calcium-52 Beam
  G4int z = 20;
  G4int a = 52;
  G4double excitEnergy = 0*MeV;
  G4double ionCharge   = 20.*eplus;
  G4double nCharge     = 0.*eplus;

  G4ParticleDefinition* particle = G4IonTable::GetIonTable()->GetIon(z,a,excitEnergy);
  G4ParticleDefinition* ntrn = G4ParticleTable::GetParticleTable()->FindParticle("neutron");

  G4double envSizeXY = 0;
  G4double envSizeZ = 0;
  G4double size = 0;
  G4double x0 = size * envSizeXY * (G4UniformRand()-0.5);
  G4double y0 = size * envSizeXY * (G4UniformRand()-0.5);
  G4double z0 = -0.5 * envSizeZ;

  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));

  for (G4int i = 0; i < events; i++){
    G4int n_particle = 1;
    fParticleGun  = new G4ParticleGun(n_particle);
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleCharge(ionCharge);
    fParticleGun->SetParticleEnergy((fragment_energy_92[i]-m_fragment) *MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(fragment_px_92[i],fragment_py_92[i],fragment_pz_92[i]));
    fParticleGun->GeneratePrimaryVertex(anEvent);
    fParticleGun->SetParticleDefinition(ntrn);
    fParticleGun->SetParticleCharge(nCharge);
    fParticleGun->SetParticleEnergy((neutron_energy_92[i]-m_neutron) *MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(neutron_px_92[i],neutron_py_92[i],neutron_pz_92[i]));
    fParticleGun->GeneratePrimaryVertex(anEvent);
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleCharge(ionCharge);
    fParticleGun->SetParticleEnergy((fragment_energy_52[i]-m_fragment) *MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(fragment_px_52[i],fragment_py_52[i],fragment_pz_52[i]));
    fParticleGun->GeneratePrimaryVertex(anEvent);
    fParticleGun->SetParticleDefinition(ntrn);
    fParticleGun->SetParticleCharge(nCharge);
    fParticleGun->SetParticleEnergy((neutron_energy_52[i]-m_neutron) *MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(neutron_px_52[i],neutron_py_52[i],neutron_pz_52[i]));
    fParticleGun->GeneratePrimaryVertex(anEvent);
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleCharge(ionCharge);
    fParticleGun->SetParticleEnergy((fragment_energy_32[i]-m_fragment) *MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(fragment_px_32[i],fragment_py_32[i],fragment_pz_32[i]));
    fParticleGun->GeneratePrimaryVertex(anEvent);
    fParticleGun->SetParticleDefinition(ntrn);
    fParticleGun->SetParticleCharge(nCharge);
    fParticleGun->SetParticleEnergy((neutron_energy_32[i]-m_neutron) *MeV);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(neutron_px_32[i],neutron_py_32[i],neutron_pz_32[i]));
    fParticleGun->GeneratePrimaryVertex(anEvent);
    //G4cout << particle->GetPDGMass() << G4endl;
  }
  //for (G4int i = 0; i < events; i++){
  //  G4int n_particle = 1;
  //  fParticleGun  = new G4ParticleGun(n_particle);
  //  fParticleGun->SetParticleDefinition(particle);
  //  fParticleGun->SetParticleCharge(ionCharge);
  //  fParticleGun->SetParticleEnergy((fragment_energy_52[i]-m_fragment) *MeV);
  //  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(fragment_px_52[i],fragment_py_52[i],fragment_pz_52[i]));
  //  fParticleGun->GeneratePrimaryVertex(anEvent);
  //  fParticleGun->SetParticleDefinition(ntrn);
  //  fParticleGun->SetParticleCharge(nCharge);
  //  fParticleGun->SetParticleEnergy((neutron_energy_52[i]-m_neutron) *MeV);
  //  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(neutron_px_52[i],neutron_py_52[i],neutron_pz_52[i]));
  //  fParticleGun->GeneratePrimaryVertex(anEvent);
  //  //G4cout << particle->GetPDGMass() << G4endl;
  //}
  //for (G4int i = 0; i < events; i++){
  //  G4int n_particle = 1;
  //  fParticleGun  = new G4ParticleGun(n_particle);
  //  fParticleGun->SetParticleDefinition(particle);
  //  fParticleGun->SetParticleCharge(ionCharge);
  //  fParticleGun->SetParticleEnergy((fragment_energy_32[i]-m_fragment) *MeV);
  //  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(fragment_px_32[i],fragment_py_32[i],fragment_pz_32[i]));
  //  fParticleGun->GeneratePrimaryVertex(anEvent);
  //  fParticleGun->SetParticleDefinition(ntrn);
  //  fParticleGun->SetParticleCharge(nCharge);
  //  fParticleGun->SetParticleEnergy((neutron_energy_32[i]-m_neutron) *MeV);
  //  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(neutron_px_32[i],neutron_py_32[i],neutron_pz_32[i]));
  //  fParticleGun->GeneratePrimaryVertex(anEvent);
  //  //G4cout << particle->GetPDGMass() << G4endl;
  //}
  vector<double> edecays_92 = edecay(events,neutrons_92,fragments_92);
  vector<double> edecays_52 = edecay(events,neutrons_52,fragments_52);
  vector<double> edecays_32 = edecay(events,neutrons_32,fragments_32);

  for (auto& edecay_92_vals : edecays_92){
    G4AnalysisManager::Instance()->FillH1(0,edecay_92_vals);
  }
  for (auto& edecay_52_vals : edecays_52){
    G4AnalysisManager::Instance()->FillH1(18,edecay_52_vals);
  }
  for (auto& edecay_32_vals : edecays_32){
    G4AnalysisManager::Instance()->FillH1(19,edecay_32_vals);
  }
  for (auto& neutron_etot_vals : neutron_energy_92){
    G4AnalysisManager::Instance()->FillH1(2,neutron_etot_vals);
  }
  for (auto& fragment_etot_vals : fragment_energy_92){
    G4AnalysisManager::Instance()->FillH1(4,fragment_etot_vals);
  }
  for (auto& neutron_px_92_vals : neutron_px_92){
    G4AnalysisManager::Instance()->FillH1(6,neutron_px_92_vals);
  }
  for (auto& neutron_py_92_vals : neutron_py_92){
    G4AnalysisManager::Instance()->FillH1(8,neutron_py_92_vals);
  }
  for (auto& neutron_pz_92_vals : neutron_pz_92){
    G4AnalysisManager::Instance()->FillH1(10,neutron_pz_92_vals);
  }
  for (auto& fragment_px_92_vals : fragment_px_92){
    G4AnalysisManager::Instance()->FillH1(12,fragment_px_92_vals);
  }
  for (auto& fragment_py_92_vals : fragment_py_92){
    G4AnalysisManager::Instance()->FillH1(14,fragment_py_92_vals);
  }
  for (auto& fragment_pz_92_vals : fragment_pz_92){
    G4AnalysisManager::Instance()->FillH1(16,fragment_pz_92_vals);
  }
  //G4cout<<G4endl;
  //G4cout<< "!@#$!@#$!@#$!@#$!@#$!@#$!#$!@#$!#$!#" << G4endl;
  //G4cout<< "This is the end" << G4endl;
  //G4cout<< "!@#$!@#$!@#$!@#$!@#$!@#$!#$!@#$!#$!#" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
