void macro_root(){
    // Read in the file that you want
    TFile* f = new TFile("Calcium.root");
    TCanvas*c1 = new TCanvas();
    // Get the relavent histograms 
    TH1D *hCalc_Ed = (TH1D*)f->Get("Caculator Edecay;1      Caculator Edecay");
    TH1D *hG4_Ed   = (TH1D*)f->Get("Geant Edecay;1  Geant Edecay");
    // Histogram characteristics
    hCalc_Ed->SetLineColor(2); // red
    hG4_Ed  ->SetLineColor(4); // blue
    hCalc_Ed->Draw();
    hG4_Ed  ->Draw("same");
    // Make a Legend
    auto legend1 = new TLegend(0.6,0.8,0.78,0.9);
    //legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    legend1->AddEntry(hCalc_Ed,"Calculator");
    legend1->AddEntry(hG4_Ed,"Geant4");
    legend1->Draw();
    ////////////////////////////
    /// Total Neutron Energy ///
    ////////////////////////////
    TCanvas*c2 = new TCanvas();
    TH1D *hCalc_Etot_n   = (TH1D*)f->Get("Caculator_Neutron_Total_Energy;1        Geant_Neutron_Total_Energy");
    TH1D *hG4_Etot_n   = (TH1D*)f->Get("Geant_Neutron_Total_Energy;1    Caculator Edecay");
    hCalc_Etot_n->SetAxisRange(1000,1200, "X");
    hG4_Etot_n->SetAxisRange(1000,1200, "X");
    hCalc_Etot_n->SetLineColor(2); // red
    hG4_Etot_n  ->SetLineColor(4); // blue
    hCalc_Etot_n->Draw();
    hG4_Etot_n  ->Draw("same");
    auto legend2 = new TLegend(0.6,0.8,0.78,0.9);
    legend2->AddEntry(hCalc_Etot_n,"Calculator");
    legend2->AddEntry(hG4_Etot_n,"Geant4");
    legend2->Draw();
    /////////////////////////////
    /// Neutron Momentum - x ///
    /////////////////////////////
    TCanvas*c3 = new TCanvas();
    TH1D *hCalc_neutron_px   = (TH1D*)f->Get("Caculator_Neutron_Px;1  Caculator_Neutron_Py");
    TH1D *hG4_px   = (TH1D*)f->Get("Geant_Neutron_Px;1      Geant_Neutron_Py");
    hCalc_neutron_px->SetAxisRange(0,100, "X");
    hG4_px->SetAxisRange(0,100, "X");
    hCalc_neutron_px->SetLineColor(2); // red
    hG4_px  ->SetLineColor(4); // blue
    hCalc_neutron_px->Draw();
    hG4_px  ->Draw("same");
    auto legend3 = new TLegend(0.6,0.8,0.78,0.9);
    legend3->AddEntry(hCalc_neutron_px,"Calculator");
    legend3->AddEntry(hG4_px,"Geant4");
    legend3->Draw();
    /////////////////////////////
    /// Neutron Momentum - y ///
    /////////////////////////////
    TCanvas*c4 = new TCanvas();
    TH1D *hCalc_neutron_py   = (TH1D*)f->Get("Caculator_Neutron_Py;1  Caculator_Neutron_Pz");
    TH1D *hG4_py   = (TH1D*)f->Get("Geant_Neutron_Py;1      Geant_Neutron_Pz");
    hCalc_neutron_py->SetAxisRange(0,100, "X");
    hG4_py->SetAxisRange(0,100, "X");
    hCalc_neutron_py->SetLineColor(2); // red
    hG4_py  ->SetLineColor(4); // blue
    hCalc_neutron_py->Draw();
    hG4_py  ->Draw("same");
    auto legend4 = new TLegend(0.6,0.8,0.78,0.9);
    legend4->AddEntry(hCalc_neutron_py,"Calculator");
    legend4->AddEntry(hG4_py,"Geant4");
    legend4->Draw();
    /////////////////////////////
    /// Neutron Momentum - z ///
    /////////////////////////////
    TCanvas*c5 = new TCanvas();
    TH1D *hCalc_neutron_pz   = (TH1D*)f->Get("Caculator_Neutron_Pz;1  Caculator Edecay");
    TH1D *hG4_pz   = (TH1D*)f->Get("Geant_Neutron_Pz;1      Geant Edecay");
    hCalc_neutron_pz->SetAxisRange(400,800, "X");
    hG4_pz->SetAxisRange(400,800, "X");
    hCalc_neutron_pz->SetLineColor(2); // red
    hG4_pz  ->SetLineColor(4); // blue
    hCalc_neutron_pz->Draw();
    hG4_pz  ->Draw("same");
    auto legend5 = new TLegend(0.6,0.8,0.78,0.9);
    legend5->AddEntry(hCalc_neutron_pz,"Calculator");
    legend5->AddEntry(hG4_pz,"Geant4");
    legend5->Draw();
    /////////////////////////////
    /// Total Fragment Energy ///
    /////////////////////////////
    TCanvas*c6 = new TCanvas();
    TH1D *hCalc_Etot_f   = (TH1D*)f->Get("Caculator_Fragment_Total_Energy;1       Geant_Fragment_Total_Energy");
    TH1D *hG4_Etot_f   = (TH1D*)f->Get("Geant_Fragment_Total_Energy;1   Caculator Edecay");
    hCalc_Etot_f->SetAxisRange(57650,57850, "X");
    hG4_Etot_f->SetAxisRange(57650,57850, "X");
    hCalc_Etot_f->SetLineColor(2); // red
    hG4_Etot_f  ->SetLineColor(4); // blue
    hCalc_Etot_f->Draw();
    hG4_Etot_f  ->Draw("same");
    auto legend6 = new TLegend(0.6,0.8,0.78,0.9);
    legend6->AddEntry(hCalc_Etot_f,"Calculator");
    legend6->AddEntry(hG4_Etot_f,"Geant4");
    legend6->Draw();
    /////////////////////////////
    /// Fragment Momentum - x ///
    /////////////////////////////
    TCanvas*c7 = new TCanvas();
    TH1D *hCalc_fragment_px   = (TH1D*)f->Get("Caculator_Fragment_Px;1 Caculator_Fragment_Py");
    TH1D *hG4_fragment_px   = (TH1D*)f->Get("Geant_Fragment_Px;1     Geant_Fragment_Py");
    hCalc_fragment_px->SetAxisRange(0,100, "X");
    hG4_fragment_px->SetAxisRange(0,100, "X");
    hCalc_fragment_px->SetLineColor(2); // red
    hG4_fragment_px  ->SetLineColor(4); // blue
    hCalc_fragment_px->Draw();
    hG4_fragment_px  ->Draw("same");
    auto legend7 = new TLegend(0.6,0.8,0.78,0.9);
    legend7->AddEntry(hCalc_fragment_px,"Calculator");
    legend7->AddEntry(hG4_fragment_px,"Geant4");
    legend7->Draw();
    /////////////////////////////
    /// Fragment Momentum - y ///
    /////////////////////////////
    TCanvas*c8 = new TCanvas();
    TH1D *hCalc_fragment_py   = (TH1D*)f->Get("Caculator_Fragment_Py;1 Caculator_Fragment_Pz");
    TH1D *hG4_fragment_py   = (TH1D*)f->Get("Geant_Fragment_Py;1     Geant_Fragment_Pz");
    hCalc_fragment_py->SetAxisRange(0,100, "X");
    hG4_fragment_py->SetAxisRange(0,100, "X");
    hCalc_fragment_py->SetLineColor(2); // red
    hG4_fragment_py  ->SetLineColor(4); // blue
    hCalc_fragment_py->Draw();
    hG4_fragment_py  ->Draw("same");
    auto legend8 = new TLegend(0.6,0.8,0.78,0.9);
    legend8->AddEntry(hCalc_fragment_py,"Calculator");
    legend8->AddEntry(hG4_fragment_py,"Geant4");
    legend8->Draw();
    /////////////////////////////
    /// Fragment Momentum - z ///
    /////////////////////////////
    TCanvas*c9 = new TCanvas();
    TH1D *hCalc_fragment_pz   = (TH1D*)f->Get("Caculator_Fragment_Pz;1 Caculator Edecay");
    TH1D *hG4_fragment_pz   = (TH1D*)f->Get("Geant_Fragment_Pz;1     Geant Edecay");
    hCalc_fragment_pz->SetAxisRange(31300,31700, "X");
    hG4_fragment_pz->SetAxisRange(31300,31700, "X");
    hCalc_fragment_pz->SetLineColor(2); // red
    hG4_fragment_pz  ->SetLineColor(4); // blue
    hCalc_fragment_pz->Draw();
    hG4_fragment_pz  ->Draw("same");
    auto legend9 = new TLegend(0.6,0.8,0.78,0.9);
    legend9->AddEntry(hCalc_fragment_pz,"Calculator");
    legend9->AddEntry(hG4_fragment_pz,"Geant4");
    legend9->Draw();
}