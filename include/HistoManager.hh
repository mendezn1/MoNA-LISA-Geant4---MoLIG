#ifndef HistoManager_h
#define HistoManager_h 1

#include "globals.hh"
#include <map>

#include "g4root.hh"
//#include "g4xml.hh"

const G4int kMaxHisto1 = 2;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class HistoManager
{
  public:
   HistoManager();
  ~HistoManager();

  private:
    void Book();

private:
    G4String fFileName;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
