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
    double intrinsicTimeResolution; // ns
    double timeStep; // ns
    bool isTimeDiscrete = false;


    //cubes
    std::map<unsigned int, std::vector<int>> cubes; // map cube number and position(first: cube number, second: x,y,z)
    std::map<unsigned int, std::vector<unsigned int>> mapCubeToChannels; // map cube number and channel number(s)

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

    TTree* configTree;
    int n_x, n_y, n_z;
    int nChannels;
    int nCubes;
    std::vector<double> *writableTimeOffsets;
    void initializeConfigTree();

    //channels
    std::map<unsigned int, std::vector<int>> channels; // map channel number and fiber position
    std::map<unsigned int, projection> channels_projection; // map channel number and fiber position
    std::map<unsigned int, double> time_offsets; // map channel number and time offset

    std::pair< unsigned int , std::vector<unsigned int> > getRandomCube();
    int getDistance(unsigned int cube, unsigned int channel);

    void setRandomTimeOffsets(double pitch);
    void setGaussianTimeOffsets(double sigma);
    void setSawToothTimeOffsets(double pitch);
    void setAlternatingTimeOffsets(double distance, double sigma);
    void resetTimeOffsets();
    void setTimeOffsets(std::string offset_type);
    void generateMapCubeToChannels();
    int verbose;

public:

    FakeSFG(int n_x, int n_y, int n_z, std::string offset_type, double timeReso, double timeStep);
    ~FakeSFG() = default;

    void printCubes();
    void printChannels();

    std::vector<FakeSFG::matchingHitsPair> simulateTimeHit();

    void generateMatchingHitsTree(int nCubeHits);
    TTree* getMatchingHitsTree();

    void setIntrinsicTimeResolution(double res);
    void setTimeStep(double ts);

    std::vector<double> getTimeOffsets(){return *writableTimeOffsets;}
    TTree* getConfigTree(){return configTree;}

    void setVerbose(int v) {verbose = v;}


};


#endif //FAKESFG_FAKESFG_HXX
