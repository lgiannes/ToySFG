#include "src/FakeSFG.hxx"


int main(int argc, char const *argv[]) {
  int n_x, n_y, n_z;
  int N;
  double timeResolution;
  std::string outputFileName;
  if(argc == 7) {
    outputFileName = argv[1];
    N = atoi(argv[2]);
    n_x = atoi(argv[3]);
    n_y = atoi(argv[4]);
    n_z = atoi(argv[5]);
    timeResolution = atof(argv[6]);
  } else{
    std::cerr<<"Usage: "<<argv[0]<<" <outputFileName> <Nhits> <n_x> <n_y> <n_z> <timeReso>"<<std::endl;
    return 1;
  }


  std::string offset_type = "alternating"; // sawtooth, gaussian, random, alternating
  FakeSFG fsfg(n_x, n_y, n_z, offset_type,timeResolution,-1);

//  fsfg.printCubes();
//  fsfg.printChannels();
  fsfg.setVerbose(0);


  fsfg.generateMatchingHitsTree(N);
  TTree* tree = fsfg.getMatchingHitsTree();
  TTree* configTree = fsfg.getConfigTree();
  auto* file = new TFile(outputFileName.c_str(), "RECREATE");
  // write n_x, n_y, n_z to the file
  configTree->Write();
  tree->Write();
  file->Write();
  file->Close();


  return 0;
}
