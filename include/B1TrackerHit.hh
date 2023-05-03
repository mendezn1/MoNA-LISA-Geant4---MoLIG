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
/// \file B1TrackerHit.hh
/// \brief Definition of the B1TrackerHit class

#ifndef B1TrackerHit_h
#define B1TrackerHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

/// Tracker hit class
///
/// It defines data members to store the trackID, chamberNb, energy deposit,
/// and position of charged particles in a selected volume:
/// - fTrackID, fChamberNB, fEdep, fPos

class B1TrackerHit : public G4VHit
{
  public:
    B1TrackerHit();
    B1TrackerHit(const B1TrackerHit&);
    virtual ~B1TrackerHit();

    // operators
    const B1TrackerHit& operator=(const B1TrackerHit&);
    G4bool operator==(const B1TrackerHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw();
    virtual void Print();

    // Set methods
    void SetTrackID  (G4int track)         { fTrackID = track; };
    void SetTrackName(G4String name)       { fParticleName = name; };
    void SetChamberNb(G4int chamb)         { fChamberNb = chamb; };
    void SetEdep     (G4double de)         { fEdep = de; };
    void SetPos      (G4ThreeVector xyz)   { fPos = xyz; };
    void SetMom      (G4ThreeVector pxpypz){ fMom = pxpypz; };

    // Get methods
    G4int GetTrackID() const     { return fTrackID; };
    G4int GetChamberNb() const   { return fChamberNb; };
    G4double GetEdep() const     { return fEdep; };
    G4ThreeVector GetPos() const { return fPos; };

  private:

      G4int         fTrackID;
      G4String      fParticleName;
      G4int         fChamberNb;
      G4double      fEdep;
      G4ThreeVector fPos;
      G4ThreeVector fMom;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<B1TrackerHit> B1TrackerHitsCollection;

extern G4ThreadLocal G4Allocator<B1TrackerHit>* B1TrackerHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* B1TrackerHit::operator new(size_t)
{
  if(!B1TrackerHitAllocator)
      B1TrackerHitAllocator = new G4Allocator<B1TrackerHit>;
  return (void *) B1TrackerHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void B1TrackerHit::operator delete(void *hit)
{
  B1TrackerHitAllocator->FreeSingle((B1TrackerHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
