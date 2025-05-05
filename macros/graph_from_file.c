


void graph_from_file(){
  TGraph * g = new TGraph("hits_stdevINOUT.txt");
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.5);
  g->SetMarkerColor(1);
  g->SetLineColor(1);
  g->Draw("APL");

  TF1 * f = new TF1("f", "[0]/sqrt(x)", 0, 1);
  f->SetParameter(0, 0.1);
  g->Fit("f");
}