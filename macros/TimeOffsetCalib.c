

void TimeOffsetCalib(std::string filename = "../output/test_Reso1000_10K.root") {
  TFile * file = new
  TFile(filename.c_str(), "READ");
  TTree *tree = (TTree *) file->Get("matchingHitsTree");

  unsigned int channel1;
  unsigned int channel2;
  double time1;
  double time2;
  double correctedTime1;
  double correctedTime2;
  double distance1;
  double distance2;
  int projection1;
  int projection2;

  tree->SetBranchAddress("channelId1", &channel1);
  tree->SetBranchAddress("channelId2", &channel2);
  tree->SetBranchAddress("time1", &time1);
  tree->SetBranchAddress("time2", &time2);
  tree->SetBranchAddress("correctedTime1", &correctedTime1);
  tree->SetBranchAddress("correctedTime2", &correctedTime2);
  tree->SetBranchAddress("distance1", &distance1);
  tree->SetBranchAddress("distance2", &distance2);
  tree->SetBranchAddress("projection1", &projection1);
  tree->SetBranchAddress("projection2", &projection2);

  int nEntries = tree->GetEntries();

  double v = 16; // cm/ns
  // distance is in cm, time is in ns

  // Open config tree
  TTree *configTree = (TTree *) file->Get("configTree");
  vector<double> *writableTimeOffsets = 0;
  int n_x, n_y, n_z;
  double timeResolution;
  double timeStep;
  configTree->SetBranchAddress("n_x", &n_x);
  configTree->SetBranchAddress("n_y", &n_y);
  configTree->SetBranchAddress("n_z", &n_z);
  configTree->SetBranchAddress("timeResolution", &timeResolution);
  configTree->SetBranchAddress("timeStep", &timeStep);
  configTree->SetBranchAddress("timeOffsets", &writableTimeOffsets);
  configTree->GetEntry(0);
  map<int, double> inputTimeOffsets;
  for (int i = 0; i < writableTimeOffsets->size(); i++) {
    inputTimeOffsets[i] = writableTimeOffsets->at(i);
//    cout << "Channel " << i << ": " << inputTimeOffsets[i] << endl;
  }
  int nChannels = inputTimeOffsets.size();

  // Start time offset calibration
  std::map<unsigned int, double> ch_sumDelta; // sum of delta
  std::map<unsigned int, double> ch_squaredDelta; // squared
  std::map<unsigned int, int> ch_entries; // entries
  std::map<unsigned int, int> ch_entries_persistent; // entries. to save after the iteration loop
  std::map<unsigned int, double> ch_offset; // sum of delta/N
  std::map<unsigned int, double> ch_offsetError; // squared/N
  std::map<unsigned int, double> lookUpTable; // final lookup table
  std::map<unsigned int, double> single_ch_offsetError; // error of each channel

  int maxIteration = 6;
  int iter = 0;
  double alpha = 0.5;
  int usedChannels = 0;

  cout<<"Start time offset calibration"<<endl;
  cout<<"Matching hits in dataset: "<<nEntries<<endl;
  cout<<"Toy SFG channels: "<<nChannels<<endl;

  while (iter < maxIteration){
    for (int hit_i = 0; hit_i<nEntries; hit_i++){
      tree->GetEntry(hit_i);

      time1 -= lookUpTable[channel1];
      time2 -= lookUpTable[channel2];

      double delta = time1 - time2 - (distance1 - distance2) / v;
      delta *= alpha;

      ch_sumDelta[channel1] += delta;
      ch_squaredDelta[channel1] += delta * delta;
      ch_entries[channel1]++;
      ch_sumDelta[channel2] -= delta;
      ch_squaredDelta[channel2] += delta * delta;
      ch_entries[channel2]++;
    } // end of loop over matching hits
    usedChannels = ch_entries.size();
    cout << "Iteration " << iter << ", used " << usedChannels << " channels." << endl;

    for (auto& entry : ch_sumDelta) {
      entry.second /= ch_entries[entry.first];
      ch_offset[entry.first] = entry.second;
      // Standard deviation is sqrt(E[X^2]-E[X]^2)
      ch_offsetError[entry.first] = sqrt( ch_squaredDelta[entry.first]/ch_entries[entry.first] - entry.second*entry.second); // std dev
      //std::cout << "Nentries is: " << ch_entries[entry.first] << std::endl;
    }

    // Fill lookup table
    for (auto& entry : ch_offset) {
//      std::cout << "Entries of channel " << entry.first << " is: " << ch_entries[entry.first] << std::endl;
      lookUpTable[entry.first] += entry.second;
    }

    // Clear the map before the next iteration
    if (iter != maxIteration-1){
      ch_sumDelta.clear();
      ch_squaredDelta.clear();
      ch_entries.clear();
    }

    iter++;
  } // end of iteration loop

  if (usedChannels < nChannels){
    int notCalibrated = nChannels - usedChannels;
    cout<<"Not calibrated channels: "<<notCalibrated<<endl;
    for(int i_ch = 0; i_ch < nChannels; i_ch++){
      if (ch_entries.find(i_ch) == ch_entries.end()){
        cout<<i_ch<<" ";
      }
    }
    cout<<endl;
  }

  TH1F* h_offset_difference = new TH1F("h_offset_difference", "Difference between input and output offsets", 100, -3, 3);
  // Print out lookup table and compare to the input
  for(int channel = 0; channel < nChannels; channel++){
//    cout<<"Channel "<<channel<<": "<<lookUpTable[channel]<<" input: "<<inputTimeOffsets[channel]<<""<<endl;
    h_offset_difference->Fill(lookUpTable[channel] - inputTimeOffsets[channel]);
  }

  TH1D* h_deltaBeforeCalib = new TH1D("h_deltaBeforeCalib", "Delta before calibration", 120, -4, 4);
  TH1D* h_deltaAfterCalib = new TH1D("h_deltaAfterCalib", "Delta after calibration", 120, -4, 4);
  for (int hit_i = 0; hit_i<nEntries; hit_i++) {
    tree->GetEntry(hit_i);
    double delta = time1 - time2 - (distance1 - distance2) / v;
    double correctedDelta = time1 - time2 - (distance1 - distance2) / v - (lookUpTable[channel1] - lookUpTable[channel2]);
    h_deltaBeforeCalib->Fill(delta);
    h_deltaAfterCalib->Fill(correctedDelta);
  }

  TCanvas* c_offDiff = new TCanvas("c_offDiff", "c_offDiff", 1800, 600);
  c_offDiff->Divide(2,1);
  c_offDiff->cd(1);
  h_offset_difference->Draw();
  h_offset_difference->Fit("gaus","Q");
  double sigma_offset = h_offset_difference->GetFunction("gaus")->GetParameter(2);
  ofstream ofs;
  ofs.open("Reso_hits_sigmaOffsets.txt", std::ofstream::out | std::ofstream::app);
  ofs << timeResolution << " " << nEntries << " " << sigma_offset << endl;
  c_offDiff->cd(2);
  h_deltaAfterCalib->SetLineColor(kBlue);
  h_deltaAfterCalib->Draw();
  h_deltaAfterCalib->Fit("gaus","Q");
  h_deltaBeforeCalib->SetLineColor(kRed);
  h_deltaBeforeCalib->Draw("same");
  h_deltaBeforeCalib->Fit("gaus","Q");
  double sigma_before = h_deltaBeforeCalib->GetFunction("gaus")->GetParameter(2);
  double sigma_after = h_deltaAfterCalib->GetFunction("gaus")->GetParameter(2);
  TLatex latex;
  latex.SetTextSize(0.03);
  latex.SetNDC();
  latex.DrawLatex(0.2, 0.91, Form("Before calib sigma: %.2f", sigma_before));
  latex.DrawLatex(0.2, 0.85, Form("After calib sigma: %.2f", sigma_after));
  TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->AddEntry(h_deltaBeforeCalib, "Before calibration", "l");
  leg->AddEntry(h_deltaAfterCalib, "After calibration", "l");
  leg->Draw();
  string reducedFilename = filename.substr(filename.find_last_of("/")+1);
  h_deltaAfterCalib->SetTitle(reducedFilename.c_str());
  h_deltaBeforeCalib->SetTitle(reducedFilename.c_str());

  c_offDiff->SaveAs(Form("TimeOffsetCalib_%s.png", reducedFilename.c_str()));


}