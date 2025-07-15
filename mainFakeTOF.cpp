//
// Created by Vedantha Srinivas Kasturi on 15.04.2025.
//

#include "src/FakeTOF.hxx"

int main(int argc, char const *argv[]) {
  int N, Nbars;
  double timeResolution;
  std::string offset_type;
  std::string outputFileName;
  if(argc == 6) {
    outputFileName = argv[1];
    N = atoi(argv[2]);
    Nbars = atoi(argv[3]); // Number of bars is equal to the number of hits
    timeResolution = atof(argv[4]);
    offset_type = argv[5];
  } else{
    std::cerr<<"Usage: "<<argv[0]<<" <outputFileName> <Nhits> <bars per panel> <timeReso (ps)> <Offset type>\n-Offset type: random, gaussian, sawtooth or alternating"<<std::endl;

    return 1;
  }

//    std::string offset_type = "random"; // sawtooth, gaussian, random, alternating
    FakeTOF ftof(offset_type,(timeResolution/1000), Nbars);
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

