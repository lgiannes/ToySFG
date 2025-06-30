#include "common_utils.h"

double ComputeCrossingMatrix(std::string filename = "/Users/lorenzo/cernbox/matchinghits/sfgd_cosmics_mHits_inclined.root", std::string title="") {

  TFile *file = new
  TFile(filename.c_str(), "READ");
  TTree *tree = (TTree *) file->Get("matchingHitsTree");
  if (tree == nullptr) {
    std::cerr << "Error: Tree matchingHitsTree not found in file " << filename << std::endl;
    return -1;
  }

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

  int nChannels = 100; // number of channels

  typedef vector< vector<double> > Matrix;
  typedef vector<double> Row;


  Matrix C;
  C.reserve(nChannels); // reserve space for nChannels rows
  std::cout<<"Setting up the matrix...\n";
  // progress bar

  for(size_t i = 0; i < nChannels; ++i)
  {
    Row row;
    row.reserve(nChannels); // reserve space for nChannels columns
    std::cout<<"Row "<<i<<std::endl;
    for(size_t j = 0; j < nChannels; ++j)
    {
      row[j] = 0;
    }

    C.push_back(row); // push each row after you fill it
  }
  std::cout<<"Finished setting up the matrix\n";

  std::cout<<"Filling the matrix...\n";
  for (int hit_i = 0; hit_i<10000; hit_i++){
    tree->GetEntry(hit_i);

    channel1 -= 2457862144;
    channel2 -= 2457862144;
    if(channel1>nChannels || channel2>nChannels) continue;

    C[channel1][channel2]++;
    C[channel2][channel1]++;

  }
  cout<<"Finished filling the matrix\n";

  cout<<"Normalizing the matrix...\n";
  // Normalize lines
  double sum[nChannels];
  for(int i=0;i<nChannels;i++) {
    sum[i]=0;
    for (int j = 0; j < nChannels; j++) {
      sum[i]+=C[i][j];
    }
  }
  for(int i=0;i<nChannels;i++) {
    for (int j = 0; j < nChannels; j++) {
      C[i][j]/=sum[i];
    }
  }
  cout<<"Finished normalizing the matrix\n";


  // Draw the matrix
  common_utils::draw_matrix(C, nChannels, title);

  return 0;
}
