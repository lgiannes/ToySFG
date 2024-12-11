
void readMatchingHits(std::string filename = "../build/test.root"){
   TFile* file = new TFile(filename.c_str(), "READ");
   TTree* tree = (TTree*) file->Get("matchingHitsTree");

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

    TH1D* hDelta = new TH1D("hDelta", "Delta", 100, -3, 3);

    for(int i=0;i<nEntries;i++){
      tree->GetEntry(i);
      double delta = time1 - time2 - (distance1 - distance2)/v;

      hDelta->Fill(delta);

      cout<<"--------------------------\nEntry "<<i<<":\n";
      cout<<"  Channel 1: "<<channel1<<", Channel 2: "<<channel2<<endl;
      cout<<"  Time 1: "<<time1<<", Time 2: "<<time2<<endl;
      cout<<"  Distance 1: "<<distance1<<", Distance 2: "<<distance2<<endl;
      cout<<"  Delta: "<<delta<<endl;
    }

    // Open config tree
    TTree* configTree = (TTree*) file->Get("configTree");
    vector<double> *writableTimeOffsets = 0;
    int n_x, n_y, n_z;
    configTree->SetBranchAddress("n_x", &n_x);
    configTree->SetBranchAddress("n_y", &n_y);
    configTree->SetBranchAddress("n_z", &n_z);
    configTree->SetBranchAddress("timeOffsets", &writableTimeOffsets);
    configTree->GetEntry(0);

    TH1D* hTimeOffsets = new TH1D("hTimeOffsets", "Time offsets", 100, -3, 3);

    int i = 0;
    for(auto to : *writableTimeOffsets){
      cout<<"Ch: "<<i<<" offset: "<< to<<endl;
      hTimeOffsets->Fill(to);
      i++;
    }

    TCanvas* c = new TCanvas("c", "c",1800, 600);
    c->Divide(2,1);
    c->cd(1);
    hDelta->Draw();
    c->cd(2);
    hTimeOffsets->Draw();





}