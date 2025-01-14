
double ComputeCrossingMatrix(std::string filename = "../output/test_3x2x2.root") {

  TFile *file = new
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
    cout << "Channel " << i << ": " << inputTimeOffsets[i] << endl;
  }
  int nChannels = inputTimeOffsets.size();

  typedef vector< vector<double> > Matrix;
  typedef vector<double> Row;

  Matrix C;
  for(size_t i = 0; i < nChannels; ++i)
  {
    Row row(nChannels);

    for(size_t j = 0; j < nChannels; ++j)
    {
      row[j] = 0;
    }

    C.push_back(row); // push each row after you fill it
  }

  for (int hit_i = 0; hit_i<nEntries; hit_i++){
    tree->GetEntry(hit_i);

    double delta_t_1 = time1 - distance1/v;
    double delta_t_2 = time2 - distance2/v;

    C[channel1][channel2]++;
    C[channel2][channel1]++;

  }

  for(int i=0;i<nChannels;i++) {
    cout<<"\n";
    for (int j = 0; j < nChannels; j++) {
      cout<<std::fixed<<std::setprecision(0)<<C[i][j]<<" ";
    }
  }
  cout<<"\n";


  // Normalize lines
  double sum[nChannels];
  for(int i=0;i<nChannels;i++) {
    sum[i]=0;
    for (int j = 0; j < nChannels; j++) {
      sum[i]+=C[i][j];
    }
  }
  cout<<"C:";
  for(int i=0;i<nChannels;i++) {
    cout<<"\n";
    for (int j = 0; j < nChannels; j++) {
      C[i][j]/=sum[i];
      cout<<std::fixed<<std::setprecision(3)<<C[i][j]<<" ";
    }
  }
  cout<<"\n";

  double a = 0.5;

  Matrix M;
  for(size_t i = 0; i < nChannels; ++i)
  {
    Row row(nChannels);

    for(size_t j = 0; j < nChannels; ++j)
    {
      if(i==j) row[j] = 1-a;
      else row[j] = a*C[i][j];
    }

    M.push_back(row); // push each row after you fill it
  }

  cout<<"M: ";
  for(int i=0;i<nChannels;i++) {
    cout<<"\n";
    for (int j = 0; j < nChannels; j++) {
      cout<<std::fixed<<std::setprecision(3)<<M[i][j]<<" ";
    }
  }
  cout<<"\n";

  double Ratio = 0; // sqrt of sum of squared magnitudes of the rows of a matrix

  TMatrixDSym TC(nChannels);
  TMatrixDSym TM(nChannels);
  for(int i=0;i<nChannels;i++) {
    double M_i_squared = 0;
    for (int j = 0; j < nChannels; j++) {
      TC[i][i]=C[i][j];
      TM[i][j]=M[i][j];
      M_i_squared += M[i][j]*M[i][j];
    }
    Ratio += M_i_squared;
    cout<<Form("||M_%d||= ",i)<<sqrt(M_i_squared)<<endl;
  }
//  TM.Print();

  Ratio = sqrt(Ratio);
  double det = TM.Determinant();
  cout<<std::scientific<<"Determinant: "<<det<<endl;

  cout<<"Ratio: "<<Ratio<<endl;



  return 0;
}
