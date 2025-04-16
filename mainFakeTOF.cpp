//
// Created by Vedantha Srinivas Kasturi on 15.04.2025.
//

#include "src/FakeTOF.hxx"

int main(int argc, char const *argv[]) {
  int N;
  double timeResolution;
  std::string outputFileName;
  if(argc == 4) {
    outputFileName = argv[1];
    N = atoi(argv[2]);
    timeResolution = atof(argv[3]);
  } else{
    std::cerr<<"Usage: "<<argv[0]<<" <outputFileName> <Nhits> <timeReso (ps)>"<<std::endl;
    return 1;
  }

    std::string offset_type = "random"; // sawtooth, gaussian, random, alternating
    FakeTOF ftof(offset_type,(timeResolution/1000));
    ftof.setVerbose(false);

    ftof.generateMatchingHitsTree(N);
    TTree* tree = ftof.getMatchingHitsTree();

    TTree* configTree = ftof.getConfigTree();
    auto* file = new TFile(outputFileName.c_str(), "RECREATE");
    // write n_x, n_y, n_z to the file
    configTree->Write();
    tree->Write();
    file->Write();
    file->Close();

    return 0;
}

