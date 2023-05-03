#include "HistoManager.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::HistoManager()
  : fFileName("Calcium")
{
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HistoManager::~HistoManager()
{
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HistoManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in HistoManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetFileName(fFileName);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);     //enable inactivation of histograms

  // Define histograms start values
  const G4String id[] = {"Caculator Edecay 92","Geant Edecay","Caculator_Neutron_Total_Energy","Geant_Neutron_Total_Energy","Caculator_Fragment_Total_Energy","Geant_Fragment_Total_Energy","Caculator_Neutron_Px","Geant_Neutron_Px","Caculator_Neutron_Py","Geant_Neutron_Py","Caculator_Neutron_Pz","Geant_Neutron_Pz","Caculator_Fragment_Px","Geant_Fragment_Px","Caculator_Fragment_Py","Geant_Fragment_Py","Caculator_Fragment_Pz","Geant_Fragment_Pz","Calculator Edecay 52","Calculator Edecay 32"};

  const G4String title[] =
       { "Caculator Edecay 92",                //0
         "Geant Edecay",                    //1
         "Caculator_Neutron_Total_Energy",  //2
         "Geant_Neutron_Total_Energy",      //3
         "Caculator_Fragment_Total_Energy", //4
         "Geant_Fragment_Total_Energy",     //5
         "Caculator_Neutron_Px",            //6
         "Geant_Neutron_Px",                //7
         "Caculator_Neutron_Py",            //8
         "Geant_Neutron_Py",                //9
         "Caculator_Neutron_Pz",            //10
         "Geant_Neutron_Pz",                //11
         "Caculator_Fragment_Px",           //12
         "Geant_Fragment_Px",               //13
         "Caculator_Fragment_Py",           //14
         "Geant_Fragment_Py",               //15
         "Caculator_Fragment_Pz",           //16
         "Geant_Fragment_Pz",               //17
         "Calculator Edecay 52"             //18
         "Calculator Edecay 32"             //19
         };

  // Default values (to be reset via /analysis/h1/set command)
  G4int nbins = 100;
  G4double vmin = 0.;
  G4double vmax = 5.;

  // Create all histograms as inactivated
  // as we have not yet set nbins, vmin, vmax
  for (G4int k=0; k<kMaxHisto1; k++) {
    G4int ih = analysisManager->CreateH1(id[k], id[k], nbins, vmin, vmax);
    analysisManager->SetH1Activation(ih, true);
  }
  for (G4int k=2; k<4; k++) {
    G4int ih = analysisManager->CreateH1(id[k], id[k], 2000, 0, 2000);
    analysisManager->SetH1Activation(ih, true);
  }
  for (G4int k=4; k<6; k++) {
    G4int ih = analysisManager->CreateH1(id[k], id[k], 60000, 0, 60000);
    analysisManager->SetH1Activation(ih, true);
  }
  for (G4int k=6; k<12; k++) {
    G4int ih = analysisManager->CreateH1(id[k], id[k], 1000, 0, 1000);
    analysisManager->SetH1Activation(ih, true);
  }
  for (G4int k=12; k<18; k++) {
    G4int ih = analysisManager->CreateH1(id[k], id[k], 40000, 0, 40000);
    analysisManager->SetH1Activation(ih, true);
  }
  for (G4int k=18; k<20; k++) {
    G4int ih = analysisManager->CreateH1(id[k], id[k], 100, 0., 5.);
    analysisManager->SetH1Activation(ih, true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
