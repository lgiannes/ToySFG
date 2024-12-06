//
// Created by Lorenzo Giannessi on 26.11.2024.
//

#include <map>
#include <vector>
#include <iostream>
#include <string>

#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>

#ifndef FAKESFG_FAKESFG_HXX
#define FAKESFG_FAKESFG_HXX


class FakeSFG {
private:

    double speed_of_light = 16; // cm/ns
    double cubeWidth = 1.027; // cm
    double intrinsicTimeResolution = 0.; // ns
    double timeStep = 2.5; // ns
    bool isTimeDiscrete = false;

    //cubes
    std::map<unsigned int, std::vector<int>> cubes; // map cube number and position(first: cube number, second: x,y,z)

    TRandom3 rnd;

    enum projection {X,Y,Z};

    struct matchingHitsPair {
        unsigned int channel1;
        unsigned int channel2;
        double time1;
        double time2;
        double correctedTime1;
        double correctedTime2;
        double distance1;
        double distance2;
        projection projection1;
        projection projection2;
    };

    TTree* matchingHitsTree;
    matchingHitsPair thisMatchingHitsPair;
    void initializeTree();
    void resetTree();
    bool isTreeFilled;

    //channels
    std::map<unsigned int, std::vector<int>> channels; // map channel number and fiber position
    std::map<unsigned int, projection> channels_projection; // map channel number and fiber position

    std::pair< unsigned int , std::vector<unsigned int> > getRandomCube();
    int getDistance(unsigned int cube, unsigned int channel);

    int verbose;

public:

    FakeSFG(int n_x, int n_y, int n_z);
    ~FakeSFG() = default;

    void printCubes();
    void printChannels();

    std::vector<FakeSFG::matchingHitsPair> simulateTimeHit();

    void generateMatchingHitsTree(int nCubeHits);
    TTree* getMatchingHitsTree();

    void setIntrinsicTimeResolution(double res);
    void setTimeStep(double ts);

    void setVerbose(int v) {verbose = v;}


};


#endif //FAKESFG_FAKESFG_HXX
