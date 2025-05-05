
// Plots the ratio of the last iteration (50) to the second last iteration (49) in the case
// of Reso=600 ps, and calibration performed with 80K hits. mini SFG (5x20x20: 600 channels).
// This is the average ratio for all channels
// It doesn't really seem to depend on the number of hits, provided it's large enough.
// It doesn't seem to depend on the resolution as well
// But only on alpha

void dummygraph(){
//  alpha	ratio	(2-α)/2
//  0.1	0.92	0.95
//  0.2	0.8996	0.9
//  0.25	0.875	0.875
//  0.3	0.8515	0.85
//  0.4	0.8052	0.8
//  0.5	0.757	0.75
//  0.75	0.643	0.625
//  0.8	0.6196	0.6
//  0.9	0.574	0.55

  TGraph * g = new TGraph();
//  g->SetPoint(0, 0.1, 0.92);
  g->SetPoint(g->GetN(), 0.2, 0.8996);
  g->SetPoint(g->GetN(), 0.25, 0.875);
  g->SetPoint(g->GetN(), 0.3, 0.8515);
  g->SetPoint(g->GetN(), 0.4, 0.8052);
  g->SetPoint(g->GetN(), 0.5, 0.757);
  g->SetPoint(g->GetN(), 0.75, 0.643);
  g->SetPoint(g->GetN(), 0.8, 0.6196);
  g->SetPoint(g->GetN(), 0.9, 0.574);

  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.5);
  g->SetMarkerColor(1);
  g->SetLineColor(1);
  g->SetTitle("Ratio of the last iteration to the second last iteration;#alpha;ratio");
  g->Draw("APL");
  g->Fit("pol1");

  TF1 * f = new TF1("f", "(2-x)/2", 0, 1);
  f->SetLineColor(kBlue);
  f->Draw("same");
  TLegend * leg = new TLegend(0.1,0.7,0.48,0.9);
  leg->AddEntry(g, "Data", "l");
  leg->AddEntry(f, "(2-#alpha)/2", "l");
  leg->Draw();

}