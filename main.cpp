#include "src/FakeSFG.hxx"


int main(int argc, char const *argv[]) {
  int n_x, n_y, n_z;
  int N;
  std::string outputFileName;
  if(argc == 6) {
    outputFileName = argv[1];
    N = atoi(argv[2]);
    n_x = atoi(argv[3]);
    n_y = atoi(argv[4]);
    n_z = atoi(argv[5]);
  } else{
    std::cerr<<"Usage: "<<argv[0]<<" <outputFileName> <Nhits> <n_x> <n_y> <n_z>"<<std::endl;
    return 1;
  }

  FakeSFG fsfg(n_x, n_y, n_z);

  fsfg.printCubes();
  fsfg.printChannels();
  fsfg.setVerbose(0);

  fsfg.generateMatchingHitsTree(N);
  TTree* tree = fsfg.getMatchingHitsTree();
  auto* file = new TFile(outputFileName.c_str(), "RECREATE");
  tree->Write();
  file->Write();
  file->Close();


  return 0;
}
