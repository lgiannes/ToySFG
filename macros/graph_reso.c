

void graph_reso(string inputfile="Reso_hits_sigmaOffsets.txt"){

  map<double,TGraph*> graphs_sigma_vs_nhits;


  ifstream ifs(inputfile);
  string line;
  while(getline(ifs, line)){
    istringstream iss(line);
    double timeReso, nhits, sigma;
    iss >> timeReso >> nhits >> sigma;
    if(graphs_sigma_vs_nhits.find(timeReso) == graphs_sigma_vs_nhits.end()){
      graphs_sigma_vs_nhits[timeReso] = new TGraph();
    }
    graphs_sigma_vs_nhits[timeReso]->SetPoint(graphs_sigma_vs_nhits[timeReso]->GetN(), nhits, sigma);
  }


  TF1* sqrtfunction = new TF1("sqrtfunction", "sqrt([1]+[0]/x)", 30000, 300000);
  sqrtfunction->SetParameters(1, 1);


  TMultiGraph* mg = new TMultiGraph();
  TLegend* leg = new TLegend(0.1,0.7,0.48,0.9);
  int i=0;
  for(auto& graph : graphs_sigma_vs_nhits){
    if(i==4) i++; // don't want yellow >:(
    graph.second->SetMarkerStyle(20);
    graph.second->SetMarkerSize(0.5);
    graph.second->SetMarkerColor(i+1);
    graph.second->SetLineColor(i+1);
    graph.second->Fit(sqrtfunction);
    graph.second->GetFunction("sqrtfunction")->SetLineColor(i+1);
    mg->Add(graph.second);
    leg->AddEntry(graph.second, Form("#sigma_{offset} = %.2f ns", graph.first), "l");
    graph.second->Draw("APL");
    i++;
  }

  mg->SetTitle("Resolution vs. number of hits;Number of hits;#sigma_{offset} [ns]");
  mg->Draw("APL");
  leg->Draw();


}