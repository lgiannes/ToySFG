

void overlay_graphs(){

  TFile * f1 = new TFile("g_W_minus_T0_0.75.root", "READ");
  TGraph* g1 = (TGraph*) f1->Get("g_W_minus_T0");
  TFile * f2 = new TFile("g_W_minus_T0_0.5.root", "READ");
  TGraph* g2 = (TGraph*) f2->Get("g_W_minus_T0");
  TFile * f3 = new TFile("g_W_minus_T0_0.25.root", "READ");
  TGraph* g3 = (TGraph*) f3->Get("g_W_minus_T0");

  TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  g1->SetMarkerStyle(20);
  g1->SetMarkerSize(0.5);
  g1->SetMarkerColor(1);
  g1->SetLineColor(1);
  g1->SetTitle("Magnitude of residual offset vector;Iteration (k);|W^{k} - T^{0}| [ns]");
  g1->Draw("APL");

  g2->SetMarkerStyle(20);
  g2->SetMarkerSize(0.5);
  g2->SetMarkerColor(2);
  g2->SetLineColor(2);
  g2->Draw("PL");

  g3->SetMarkerStyle(20);
  g3->SetMarkerSize(0.5);
  g3->SetMarkerColor(3);
  g3->SetLineColor(3);
  g3->Draw("PL");

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(g1, "#alpha=0.75", "l");
  leg->AddEntry(g2, "#alpha=0.5", "l");
  leg->AddEntry(g3, "#alpha=0.25", "l");
  leg->Draw("same");
  gPad->SetLogy();

}