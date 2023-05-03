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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"
#include "B1TrackerSD.hh" // NM added 2023-3-19

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh" // NM added 2023-3-19



#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"

#include "G4UserLimits.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  G4int ncomponents, natoms;
  G4double a, z, density;
  G4double temperature, pressure, volume, gas_constant, mass;

  // Nitrogen (vacuum)
  a = 14.01*g/mole;
  temperature = 239.; // kelvin
  pressure = 0.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001/760; // atm
  volume = 0.012; //m3
  gas_constant = 8.2057; // m3*atm/mol*K
  mass = ((pressure*volume)/(gas_constant*temperature))*a;
  density = mass/volume;
  G4Material* elN = new G4Material("Nitrogen", z=7., a, density);

  //
  // WORLD
  //
  G4double world_sizeXY = 25*cm;
  G4double world_sizeZ  = 25*cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");

  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        elN,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  // SENSITIVE DETECTOR
  //
  G4double sensdet_sizeXY = 25.0*cm;
  G4double sensdet_sizeZ  = 1.0*cm;
  G4Material* sensdet_mat = nist->FindOrBuildMaterial("G4_Galactic");
  G4ThreeVector posSensdet = G4ThreeVector(0.*cm, 0.*cm, 12*cm);

  G4Box* sensdet =
    new G4Box("SensDet",                       //its name
              0.5*sensdet_sizeXY, 
              0.5*sensdet_sizeXY, 
              0.5*sensdet_sizeZ);     //its size

  G4LogicalVolume* logicSensdet =
    new G4LogicalVolume(sensdet,          //its solid
                        elN,           //its material
                        "SensDet");            //its name

  new G4PVPlacement(0,                     //no rotation
                      posSensdet,       //at (0,0,0)
                      logicSensdet,            //its logical volume
                      "SensDet",               //its name
                      logicWorld,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

  //
  // TARGET
  //
  //G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_Be");
  //G4ThreeVector pos2 = G4ThreeVector(0*cm, 0*cm, 1*cm);
//
  //// Target Dimmensions
  //G4double shape2_dx = 2.80*cm;
  //G4double shape2_dy = 2.80*cm;
  //G4double shape2_dz = 0.41*cm;
//
  //G4Box* solidShape2 =
  //  new G4Box("Shape2",                      //its name
  //            0.5*shape2_dx,
  //            0.5*shape2_dy,
  //            0.5*shape2_dz); //its size
//
  //G4LogicalVolume* logicShape2 =
  //  new G4LogicalVolume(solidShape2,         //its solid
  //                      shape2_mat,          //its material
  //                      "Shape2");           //its name
//
  //new G4PVPlacement(0,                       //no rotation
  //                  pos2,                    //at position
  //                  logicShape2,             //its logical volume
  //                  "Shape2",                //its name
  //                  logicWorld,              //its mother  volume
  //                  false,                   //no boolean operation
  //                  0,                       //copy number
  //                  checkOverlaps);          //overlaps checking
//
  // Set Shape2 as scoring volume
  //
  fScoringVolume = logicSensdet;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1DetectorConstruction::ConstructSDandField()
{
  // Sensitive detectors

  G4String trackerChamberSDname = "B1/TrackerChamberSD";
  B1TrackerSD* aTrackerSD = new B1TrackerSD(trackerChamberSDname,
                                            "TrackerHitsCollection");
  G4SDManager::GetSDMpointer()->AddNewDetector(aTrackerSD);
  // Setting aTrackerSD to all logical volumes with the same name 
  // of "SensDet".
  SetSensitiveDetector("SensDet", aTrackerSD, true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......