//
// Created by Vedantha Srinivas Kasturi on 14.04.2025.
//

#ifndef FAKESFG_FAKETOF_HXX
#define FAKESFG_FAKETOF_HXX

#include <map>
#include <vector>
#include <iostream>
#include <string>

#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include "TVector3.h"


class FakeTOF {

private:
    double speed_of_light = 30; // cm/ns
    double intrinsicTimeResolution; // ns
    double barLength = 200; // cm
    double barSeparation = 15; // cm

    int barsPerPanel = 10; // Number of bars per panel

    double TopBar0[3] = {0, 100, -75}; // From upstream to downstream
    double BottomBar0[3] = {0, -100, -75}; // From upstream to downstream
    double UpstreamBar0[3] = {0, 75, -100}; // From top to bottom
    double DownstreamBar0[3] = {0, 75, 100}; // From top to bottom
    double NorthBar0[3] = {-100, 75, 0}; // From top to bottom
    double SouthBar0[3] = {100, 75, 0}; // From top to bottom

    std::map<int, std::vector<double>> Topbars; // map bar number and center position(first: bar number, second: x,y,z)
    std::map<int, std::vector<double>> Bottombars; // map bar number and center position(first: bar number, second: x,y,z)
    std::map<int, std::vector<double>> Upstreambars; // map bar number and center position(first: bar number, second: x,y,z)
    std::map<int, std::vector<double>> Downstreambars; // map bar number and center position(first: bar number, second: x,y,z)
    std::map<int, std::vector<double>> Northbars; // map bar number and center position(first: bar number, second: x,y,z)
    std::map<int, std::vector<double>> Southbars; // map bar number and center position(first: bar number, second: x,y,z)

    std::map<int, double> time_offsetsTop; // map bar number and time offset
    std::map<int, double> time_offsetsBottom; // map bar number and time offset
    std::map<int, double> time_offsetsUpstream; // map bar number and time offset
    std::map<int, double> time_offsetsDownstream; // map bar number and time offset
    std::map<int, double> time_offsetsNorth; // map bar number and time offset
    std::map<int, double> time_offsetsSouth; // map bar number and time offset

    // Enum for panel numbers

    enum PanelId {
        TOP = 0,
        BOTTOM = 1,
        UPSTREAM = 2,
        DOWNSTREAM = 3,
        NORTH = 4,
        SOUTH = 5
    };




    TRandom3 rnd;

    struct matchingHitsPair {
        int panel1;
        int panel2;
        int bar1;
        int bar2;
        double distance;
        double flightTimePostOffset;
        double flightTimeRaw;
        double time1;
        double time2;
        double correctedTime1;
        double correctedTime2;
        double hit1Position[3];
        double hit2Position[3];
        double DeltaTimePostOffset;
    };

    TTree* matchingHitsTree;
    matchingHitsPair thisMatchingHitsPair;

    TTree* configTree;

    std::vector<double> *writableTimeOffsetsTop;
    std::vector<double> *writableTimeOffsetsBottom;
    std::vector<double> *writableTimeOffsetsUpstream;
    std::vector<double> *writableTimeOffsetsDownstream;
    std::vector<double> *writableTimeOffsetsNorth;
    std::vector<double> *writableTimeOffsetsSouth;

    void initializeTree();
    void initializeConfigTree();

    void setTimeOffsets(std::string offset_type);

    void setRandomTimeOffsets(double pitch);
    void setGaussianTimeOffsets(double sigma);
    void setSawToothTimeOffsets(double pitch);
    void setAlternatingTimeOffsets(double distance, double sigma);
    void resetTimeOffsets();
    bool verbose;



public:
    FakeTOF(std::string offset_type, double timeReso); // Initialise the ToF detector
    ~FakeTOF() = default;

    FakeTOF::matchingHitsPair simulateTimeHit(); // Simulate a time hit

    void generateMatchingHitsTree(int NHits); // Generate a tree with matching hits

    void setVerbose(bool v) { verbose = v; } // Set verbosity level

    TTree* getMatchingHitsTree() { return matchingHitsTree; } // Get the matching hits tree
    TTree* getConfigTree() { return configTree; } // Get the configuration tree


};


#endif //FAKESFG_FAKETOF_HXX
