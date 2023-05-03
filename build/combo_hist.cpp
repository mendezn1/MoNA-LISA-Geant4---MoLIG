void combo_hist(){
    TFile* f = new TFile("Calcium.root");
    TCanvas*c1 = new TCanvas();
    TH1D *hCalc_Ed_92 = (TH1D*)f->Get("Caculator Edecay 92;1   Caculator Edecay 92");
    TH1D *hCalc_Ed_52 = (TH1D*)f->Get("Calculator Edecay 52;1  Calculator Edecay 52");
    TH1D *hCalc_Ed_32 = (TH1D*)f->Get("Calculator Edecay 32;1  Calculator Edecay 32");
    TH1D *hGeant      = (TH1D*)f->Get("Geant Edecay;1  Geant Edecay");
    

    auto hedecay_tot = new TH1D(*hCalc_Ed_92);
    hedecay_tot->Add(hCalc_Ed_52);
    hedecay_tot->Add(hCalc_Ed_32);

    hCalc_Ed_92->SetLineColor(kRed);
    hCalc_Ed_52->SetLineColor(kGreen+2);
    hCalc_Ed_32->SetLineColor(kBlue);
    hedecay_tot->SetLineColor(kBlack);
    //hGeant     ->SetLineColor(kMagenta);

    hCalc_Ed_92->Draw("SAME H");
    hCalc_Ed_52->Draw("SAME H");
    hCalc_Ed_32->Draw("SAME H");
    hedecay_tot->Draw("SAME H");
    //hGeant     ->Draw("SAME H");

    auto legend = new TLegend(0.6,0.775,0.78,0.9);
    legend->AddEntry(hCalc_Ed_92,"9/2+");
    legend->AddEntry(hCalc_Ed_52,"5/2+");
    legend->AddEntry(hCalc_Ed_32,"3/2+");
    legend->AddEntry(hedecay_tot,"Edecay Tot");
    legend->Draw();

    // Edecay Comparison (MoLIG vs. UnNPy) //

    TCanvas*c2 = new TCanvas();

    hedecay_tot->SetLineColor(kBlue);
    hGeant->SetLineColor(kRed);

    hGeant->Draw();
    hedecay_tot->Draw("SAME");
    
    auto legend1 = new TLegend(0.6,0.775,0.78,0.9);
    legend1->AddEntry(hGeant,"Geant Edecay Tot");
    legend1->AddEntry(hedecay_tot,"UnNPy Edecay Tot");
    legend1->Draw();

    // Total Edecay Difference //

    TCanvas*c3 = new TCanvas();
    auto edecay_tot_difference = new TH1D(*hedecay_tot);
    edecay_tot_difference->Add(hGeant,-1);
    edecay_tot_difference->Draw("SAME");

    // Total Edecay Ratio //

    TCanvas*c4 = new TCanvas();
    auto edecay_tot_ratio = new TH1D(*hedecay_tot);
    edecay_tot_ratio->Divide(hGeant);
    edecay_tot_ratio->Draw("SAME");
    
}