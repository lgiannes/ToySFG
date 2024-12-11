#include "TimeOffsetCalib.c"

void sigmaOffset_vs_channels() {

  vector<int> nzChannels = {2,3,4,5,8,10,12,15};

  TGraph* g_sigma_vs_nch = new TGraph();
  TGraph* g_sigma_vs_ncubes = new TGraph();

  for (int nz : nzChannels) {
    string inputfile = Form("../output/test_%dch.root", nz);
    double nCubes = nz*(4*nz)*(4*nz);
    double nChannels = nz*(4*nz) + nz*(4*nz) + (4*nz)*(4*nz);
    double sigma = TimeOffsetCalib(inputfile,false,true);
    g_sigma_vs_nch->SetPoint(g_sigma_vs_nch->GetN(), nChannels, sigma);
    g_sigma_vs_ncubes->SetPoint(g_sigma_vs_ncubes->GetN(), nCubes, sigma);
  }

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  c1->Divide(2,1);
  c1->cd(1);
  g_sigma_vs_nch->SetMarkerStyle(20);
  g_sigma_vs_nch->SetMarkerSize(0.5);
  g_sigma_vs_nch->SetMarkerColor(1);
  g_sigma_vs_nch->SetLineColor(1);
  g_sigma_vs_nch->SetTitle("#sigma_{offset} vs. number of channels;Number of channels;#sigma_{offset} [ns]");
  g_sigma_vs_nch->Draw("APL");
  c1->cd(2);
  g_sigma_vs_ncubes->SetMarkerStyle(20);
  g_sigma_vs_ncubes->SetMarkerSize(0.5);
  g_sigma_vs_ncubes->SetMarkerColor(2);
  g_sigma_vs_ncubes->SetLineColor(2);
  g_sigma_vs_ncubes->SetTitle("#sigma_{offset} vs. number of cubes;Number of cubes;#sigma_{offset} [ns]");
  g_sigma_vs_ncubes->Draw("APL");

}

