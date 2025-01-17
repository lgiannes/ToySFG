
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

  //////////////////////////////////////// ALPHA
  double a = 0.5;
  //////////////////////////////////////////////

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


  TMatrixDSym TC(nChannels);
  TMatrixDSym TM(nChannels);
  TMatrixDSym TI(nChannels);
  TMatrixDSym TJ(nChannels); // matrix of all 1s
  for(int i=0;i<nChannels;i++) {
    double M_i_squared = 0;
    for (int j = 0; j < nChannels; j++) {
      TC[i][j]=C[i][j];
      TM[i][j]=M[i][j];
      TJ[i][j]=1;
      if(i==j) {
        TI[i][j]=1;
      }
      else {
        TI[i][j]=0;
      }
      M_i_squared += M[i][j];
    }
    cout<<Form("M_%d= ",i)<<(M_i_squared)<<"  - sum of the elements in a row of the  M matrix."<<endl;
  }
//  TM.Print();
  cout<<"I-C: ";
  for(int i=0;i<nChannels;i++) {
    cout<<"\n";
    for (int j = 0; j < nChannels; j++) {
      cout<<std::fixed<<std::setprecision(3)<<TI[i][j]-TC[i][j]<<" ";
    }
  }
  cout<<"\n";

  TMatrixD TA(TJ, TMatrixD::kMult, TI-TC); // matrix correction to implement the averaging to zero
  for(int i=0;i<nChannels;i++) {
    for (int j = 0; j < nChannels; j++) {
      TA[i][j] *= a/nChannels;
    }
  }
  cout<<"A: ";
  for(int i=0;i<nChannels;i++) {
    cout<<"\n";
    for (int j = 0; j < nChannels; j++) {
      cout<<std::fixed<<std::setprecision(3)<<TA[i][j]<<" ";
    }
  }
  cout<<"\n";

  TMatrixD TMA = TA + TM;
  TMatrixDEigen TM_eigen_corrected(TMA);
  TVectorD eigenvalues_corrected(nChannels);
  eigenvalues_corrected = TM_eigen_corrected.GetEigenValues();
  TMatrixD eigenvectors_corrected = TM_eigen_corrected.GetEigenVectors();

  cout<<"M+A: ";
  for(int i=0;i<nChannels;i++) {
    cout<<"\n";
    double sumrow = 0;
    for (int j = 0; j < nChannels; j++) {
      cout<<std::fixed<<std::setprecision(3)<<TMA[i][j]<<" ";
      sumrow += TMA[i][j];
    }
    cout<<"  sum: "<<sumrow;
  }
  cout<<"\n";

  TMatrixDEigen TM_eigen(TM);
  TVectorD eigenvalues(nChannels);
  TVectorD eigenvalues_Im(nChannels);
  eigenvalues = TM_eigen.GetEigenValues();
  eigenvalues_Im = TM_eigen.GetEigenValuesIm();
  TMatrixD eigenvectors = TM_eigen.GetEigenVectors();
//  eigenvectors.Print();

  // Compute M' = (I-C)M(I-C)^-1
  TMatrixDSym A = TI - TC;
  TMatrixDSym Ainv = A.Invert();
  TMatrixD Mprime(A, TMatrixD::kMult, TM);
  Mprime *= Ainv;
//  Mprime.Print();
  cout<<std::fixed<<std::setprecision(4);
  cout<<"Eigenvalues:\n";
  for(int i=0;i<nChannels;i++) {
    cout<<eigenvalues_corrected[i]<<"\n";
  }

  // Print first row of the eigenvector matrix
  cout<<"Eigenvectors first column:\n";
  double sum_eigen = 0;
  for(int i=0;i<nChannels;i++) {
    sum_eigen+=eigenvectors_corrected[i][0];
  }
  for(int i=0;i<nChannels;i++) {
    cout<<eigenvectors_corrected[i][0] <<" ";
  }
  cout<<"\n";
  cout<< "time offsets: "<<endl;
  for(double deltaT_0 : *writableTimeOffsets){
    cout<<deltaT_0<<" ";
  }
  cout<<endl;


  // Compute powers of M matrix
  int K = 40;
  TMatrixD M_k(TI);
  TMatrixD M_series(TI);
  TMatrixD DeltaTk(nChannels,1);
  TMatrixD DeltaTkMinusOne(nChannels,1);
  TMatrixD DeltaT0(nChannels,1);
  TMatrixD DeltaT0_series(nChannels,1);
  TMatrixD correction_vector_k(nChannels,1);
  TMatrixD correction_vector_sum_k(nChannels,1);
  for(int i=0;i<nChannels;i++) {
    DeltaT0[i][0] = (*writableTimeOffsets)[i];
    correction_vector_sum_k[i][0] = 0;
  }
  DeltaT0_series = DeltaT0;
  DeltaTkMinusOne = DeltaT0;
  for(int k=1;k<K;k++) {
    DeltaTkMinusOne = M_k*DeltaT0;
    M_k *= TMA;
    M_series += M_k;
    DeltaTk = M_k*DeltaT0;
    DeltaT0_series += DeltaTk;
    cout<<"DeltaT^"<<k<<":\n";
    double M_series_norm = 0;
    double DeltaTk_norm = 0;
    double DeltaTk_average = 0;
    double DeltaTk_series_norm = 0;
    for(int i=0;i<nChannels;i++) {
      cout<<DeltaTk[i][0]<<" ";
      DeltaTk_norm += DeltaTk[i][0]*DeltaTk[i][0];
      DeltaTk_average += DeltaTk[i][0];
      DeltaTk_series_norm += DeltaT0_series[i][0]*DeltaT0_series[i][0];
      for (int j = 0; j < nChannels; j++) {
//        cout<<std::fixed<<std::setprecision(3)<<M_k[i][j]<<" ";
        M_series_norm += M_k[i][j]*M_k[i][j];
      }
//      cout<<"\n";
    }
    cout<<"  average: "<<DeltaTk_average<<"\n";
    cout<<"Sum of DeltaT^"<<k<<":\n";
    for(int i=0;i<nChannels;i++) {
      cout<<DeltaT0_series[i][0]<<" ";
    }
    cout<<"\n";
    correction_vector_k = DeltaTk-DeltaTkMinusOne;
    cout<<"correction vector d^"<<k<<":\n";
    for(int i=0;i<nChannels;i++) {
      cout<<correction_vector_k[i][0]<<" ";
    }
    cout<<"\n";
    correction_vector_sum_k += correction_vector_k;
    cout<<"correction vector sum d^"<<k<<":\n";
    for(int i=0;i<nChannels;i++) {
      cout<<correction_vector_sum_k[i][0]<<" ";
    }
    cout<<"\n";

    DeltaTk_average /= nChannels;
    M_series_norm = sqrt(M_series_norm);
    DeltaTk_norm = sqrt(DeltaTk_norm);
    cout<<"Norm of  M_"<<k<<": "<<M_series_norm<<endl;
    cout<<"Norm of DT_"<<k<<": "<<DeltaTk_norm<<endl;
    cout<<"Norm of DTk series: "<<sqrt(DeltaTk_series_norm)<<endl;
  }

//  double det = TM.Determinant();
//  cout<<std::scientific<<"Determinant: "<<det<<endl;

  cout<< "time offsets: "<<endl;
  for(double deltaT_0 : *writableTimeOffsets){
    cout<<deltaT_0<<" ";
  }
  cout<<endl;
  cout<<" Perron-Frobenius eigenvalue: "<<abs(eigenvalues_corrected[0])<<" (largest absolute eigenvalue)"<<endl;
  cout<<" Convergence rate: "<<abs(eigenvalues_corrected[1])<<" (second largest eigenvalue)"<<endl;



  return 0;
}
