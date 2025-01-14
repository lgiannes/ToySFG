

void graph_sigmaOffset_vs_iter(string inputfile="iter_vs_sigma_offset.txt") {

  TGraph * g_sigma_vs_iter = new
  TGraph();

  string header1;
  string header2;
  ifstream ifs(inputfile);
  string line;
  // ignore first two lines
  getline(ifs, line);
  header1 = line;
  getline(ifs, line);
  header2 = line;
  while (getline(ifs, line)) {
    istringstream iss(line);
    double iter, sigma;
    iss >> iter >> sigma;
    g_sigma_vs_iter->SetPoint(g_sigma_vs_iter->GetN(), iter, sigma);
  }

  TCanvas * c1 = new
  TCanvas("c1", "c1", 800, 600);
  g_sigma_vs_iter->SetMarkerStyle(20);
  g_sigma_vs_iter->SetMarkerSize(0.5);
  g_sigma_vs_iter->SetMarkerColor(1);
  g_sigma_vs_iter->SetLineColor(1);
  g_sigma_vs_iter->SetTitle("#sigma_{offset} vs. number of iterations;Number of iterations;#sigma_{offset} [ns]");
  g_sigma_vs_iter->Draw("APL");
  g_sigma_vs_iter->SetTitle(header1.c_str());

  // Geometrical series

}