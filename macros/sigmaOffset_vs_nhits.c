#include "TimeOffsetCalib.c"

void sigmaOffset_vs_nhits(){



  vector<int> channels = {2,3,4,5,6,8,10,12,15}; //

  TMultiGraph* mg = new TMultiGraph();
  TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);

  for(int ch : channels) {
    string inputfile = Form("../output/test_%dch.root", ch);

    vector<int> nhits = {300000, 100000, 50000, 30000}; //1500000,1000000,500000,
    TGraph * g_sigma_vs_nhits = new TGraph();

    for (int n: nhits) {
      double sigma = TimeOffsetCalib(inputfile, false, true, n);
      g_sigma_vs_nhits->SetPoint(g_sigma_vs_nhits->GetN(), n, sigma);
    }
    leg->AddEntry(g_sigma_vs_nhits, Form("toy SFG:  %dx%dx%d cubes",4*ch,ch,4*ch), "l");

    mg->Add(g_sigma_vs_nhits);
  }

  TF1* sqrtfunction = new TF1("sqrtfunction", "sqrt([1]+[0]/x)", 30000, 300000);
  sqrtfunction->SetParameters(1, 1);

  for( int i = 0; i < mg->GetListOfGraphs()->GetSize(); i++ ) {
    TGraph * g = (TGraph*) mg->GetListOfGraphs()->At(i);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.5);
    g->SetMarkerColor(i+1);
    g->SetLineColor(i+1);
    g->Fit(sqrtfunction);
    g->GetFunction("sqrtfunction")->SetLineColor(i+1);
  }

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  mg->SetTitle("#sigma_{offset} vs. number of hits;Number of hits;#sigma_{offset} [ns]");
  mg->Draw("APL");
  leg->Draw();


}

