//
// Created by Lorenzo Giannessi on 26.11.2024.
//

#include "FakeSFG.hxx"


FakeSFG::FakeSFG(int n_x, int n_y, int n_z) {

    if(n_x < 1 || n_y < 1 || n_z < 1){
        std::cerr<<"Error: n_x, n_y and n_z must be greater than 0"<<std::endl;
        return;
    }

    //initialize cubes
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
            for (int k = 0; k < n_z; k++) {
                std::vector<int> pos = {i, j, k};
                cubes[i * n_y * n_z + j * n_z + k] = pos;
            }
        }
    }
    std::cout<< cubes.size() <<" cubes initialized  --  Cross check: "<<n_x*n_y*n_z<<" == "<<cubes.size()<<std::endl;
    //initialize channels
    int index = 0;
    // Face XY (fibers along z)
    for (int i = 0; i < n_x; i++) {
        for (int j = 0; j < n_y; j++) {
          std::vector<int> pos = {i, j, -1};
          channels[index] = pos;
          channels_projection[index] = Z;
          index++;
        }
    }
    // Face XZ (fibers along y)
    for (int i = 0; i < n_x; i++) {
        for (int k = 0; k < n_z; k++) {
          std::vector<int> pos = {i, -1, k};
          channels[index] = pos;
          channels_projection[index] = Y;
          index++;
        }
    }
    // Face YZ (fibers along x)
    for (int j = 0; j < n_y; j++) {
        for (int k = 0; k < n_z; k++) {
          std::vector<int> pos = {-1, j, k};
          channels[index] = pos;
          channels_projection[index] = X;
          index++;
        }
    }

    setVerbose(0); // by default verbose is off
    std::cout<< channels.size() <<" channels initialized  --  Cross check: "<<n_x*n_y + n_y*n_z + n_z*n_x<<" == "<<channels.size()<<std::endl;

    std::cout<<"Initializing tree...";
    initializeTree();
    isTreeFilled = false;
    std::cout<<" done"<<std::endl;

}

void FakeSFG::printCubes() {
    for (auto const& [key, val] : cubes) {
        std::cout << "Cube " << key << " at position: ";
        for (auto const& v : val) {
            std::cout << v << " ";
        }
        std::cout << std::endl;
    }
}

void FakeSFG::initializeTree() {
    matchingHitsTree = new TTree("matchingHitsTree", "Tree containing the matching hits");
    matchingHitsTree->Branch("time1", &thisMatchingHitsPair.time1, "time1/D");
    matchingHitsTree->Branch("time2", &thisMatchingHitsPair.time2, "time2/D");
    matchingHitsTree->Branch("correctedTime1", &thisMatchingHitsPair.correctedTime1, "correctedTime1/D");
    matchingHitsTree->Branch("correctedTime2", &thisMatchingHitsPair.correctedTime2, "correctedTime2/D");
    matchingHitsTree->Branch("distance1", &thisMatchingHitsPair.distance1, "distance1/D");
    matchingHitsTree->Branch("distance2", &thisMatchingHitsPair.distance2, "distance2/D");
    matchingHitsTree->Branch("channelId1", &thisMatchingHitsPair.channel1, "channelId1/i");
    matchingHitsTree->Branch("channelId2", &thisMatchingHitsPair.channel2, "channelId2/i");
    matchingHitsTree->Branch("projection1", &thisMatchingHitsPair.projection1, "projection1/I");
    matchingHitsTree->Branch("projection2", &thisMatchingHitsPair.projection2, "projection2/I");
}

void FakeSFG::printChannels() {

    int x_fibers = 0;
    int y_fibers = 0;
    int z_fibers = 0;
    for (auto const& [key, val] : channels) {
      std::cout << "Channel " << key << " at position: ";
      for (auto const &v: val) {
        std::cout << v << " ";
      }
      std::cout << " -- Projection: ";
      switch (channels_projection[key]) {
        case X:
          std::cout << "X";
          x_fibers++;
          break;
        case Y:
          std::cout << "Y";
          y_fibers++;
          break;
        case Z:
          std::cout << "Z";
          z_fibers++;
          break;
      }
      std::cout << std::endl;
    }
    std::cout<<"X fibers: "<<x_fibers<<"\nY fibers: "<<y_fibers<<"\nZ fibers: "<<z_fibers<<std::endl;
}

// get a random cube and returns the three channels reading the cube + the cube number
std::pair< unsigned int , std::vector<unsigned int> > FakeSFG::getRandomCube() {
    unsigned int random = rand() % cubes.size();

    std::vector<int> pos_cube = cubes[random];
//    std::cout<<"cube: "<<random<<std::endl;
//    std::cout<<"cube pos: "<<pos_cube[0]<<" "<<pos_cube[1]<<" "<<pos_cube[2]<<std::endl;
    std::vector<unsigned int> channels_reading;
    for (auto const& [key, val] : channels) {
        std::vector<int> pos_channel = val;
        if (channels_projection[key] == X && pos_channel[1] == pos_cube[1] && pos_channel[2] == pos_cube[2]) {
            channels_reading.push_back(key);
//            std::cout<<"debug: X:"<<key<<std::endl;
//            std::cout<<pos_channel[1]<<" "<<pos_cube[1]<<" - "<<pos_channel[2]<<" "<<pos_cube[2]<<std::endl;
        }
        if (channels_projection[key] == Y && pos_channel[2] == pos_cube[2] && pos_channel[0] == pos_cube[0]) {
            channels_reading.push_back(key);
//            std::cout<<"debug: Y:"<<key<<std::endl;
//            std::cout<<pos_channel[2]<<" "<<pos_cube[2]<<" - "<<pos_channel[0]<<" "<<pos_cube[0]<<std::endl;
        }
        if (channels_projection[key] == Z && pos_channel[1] == pos_cube[1] && pos_channel[0] == pos_cube[0]) {
            channels_reading.push_back(key);
//            std::cout<<"debug: Z:"<<key<<std::endl;
//            std::cout<<pos_channel[1]<<" "<<pos_cube[1]<<" - "<<pos_channel[0]<<" "<<pos_cube[0]<<std::endl;
        }
    }
    if(channels_reading.size() != 3){
        std::cerr<<"ERROR: Cube "<<random<<" is not read by 3 channels"<<std::endl;
    }

    std::pair< unsigned int , std::vector<unsigned int> > result = std::make_pair(random, channels_reading);


    return result;

}

int FakeSFG::getDistance(unsigned int cube, unsigned int channel) {
    std::vector<int> pos_cube = cubes[cube];
    std::vector<int> pos_channel = channels[channel];

    int distance = -1;

    if (channels_projection[channel] == X && pos_channel[1] == pos_cube[1] && pos_channel[2] == pos_cube[2]) {
        distance = abs(pos_channel[0] - pos_cube[0]) ;
    }else if (channels_projection[channel] == Y && pos_channel[2] == pos_cube[2] && pos_channel[0] == pos_cube[0]) {
        distance = abs(pos_channel[1] - pos_cube[1]) ;
    }else if (channels_projection[channel] == Z && pos_channel[1] == pos_cube[1] && pos_channel[0] == pos_cube[0]) {
        distance = abs(pos_channel[2] - pos_cube[2]) ;
    }else{
        std::cerr<<"Error: Channel "<<channel<<" does not read cube "<<cube<<std::endl;
    }

    return distance;
}

std::vector<FakeSFG::matchingHitsPair> FakeSFG::simulateTimeHit() {
  if(verbose) {
    std::cout << "----------------------------------------\nDEBUG: simulating a time hit.\n";
  }
  std::pair<unsigned int, std::vector<unsigned int> > cube_reading = getRandomCube();
  std::vector<unsigned int> channels_reading = cube_reading.second;
  unsigned int cube = cube_reading.first;

  double times[3];
  double time_x = getDistance(cube, channels_reading[0]) * cubeWidth / speed_of_light;
  times[0] = time_x;
  double time_y = getDistance(cube, channels_reading[1]) * cubeWidth / speed_of_light;
  times[1] = time_y;
  double time_z = getDistance(cube, channels_reading[2]) * cubeWidth / speed_of_light;
  times[2] = time_z;

  std::vector<matchingHitsPair> ThreeHitPairs;
  for(int i = 0; i < 3; i++) {
    matchingHitsPair pair;
    pair.channel1 = channels_reading[i];
    pair.time1 = times[i];
    pair.distance1 = getDistance(cube, channels_reading[i]);
    pair.channel2 = channels_reading[(i + 1) % 3];
    pair.time2 = times[(i + 1) % 3];
    pair.distance2 = getDistance(cube, channels_reading[(i + 1) % 3]);
    ThreeHitPairs.push_back(pair);
  }

  if (verbose) {
    std::cout << "Cube: " << cube << std::endl;
    std::cout << "Cube position: ";
    for (auto const &v: cubes[cube]) {
      std::cout << v << " ";
    }
    std::cout << std::endl;
    std::cout << "Channels reading: ";
    for (auto const &v: channels_reading) {
      std::cout << v << " ";
    }
    std::cout << std::endl;
    for (auto const &v: channels_reading) {
      std::cout << "Distance from cube to channel " << v << ": " << getDistance(cube, v) << std::endl;
    }
    std::cout << "Hit in cube " << cube << " at position: ";
    for (auto const &v: cubes[cube]) {
      std::cout << v << " ";
    }
    std::cout << std::endl;
    std::cout << "Time hit in channel " << channels_reading[0] << ": " << time_x << " ns\n";
    std::cout << "Time hit in channel " << channels_reading[1] << ": " << time_y << " ns\n";
    std::cout << "Time hit in channel " << channels_reading[2] << ": " << time_z << " ns\n";
  }

  // hitPairs[0]: XY, hitPairs[1]: YZ, hitPairs[2]: XZ
  return ThreeHitPairs;
}

void FakeSFG::generateMatchingHitsTree(int nCubeHits) {
  if(!matchingHitsTree) {
    std::cerr<<"Error: Tree not initialized"<<std::endl;
    return;
  }

  for(int i = 0; i < nCubeHits; i++) {
    std::vector<matchingHitsPair> hits = simulateTimeHit();
    for (auto const &hit: hits) {
      thisMatchingHitsPair = hit;
      matchingHitsTree->Fill();
    }
  }

  std::cout << "Simulated " << nCubeHits << " cube hits.\n";
  std::cout << "Output TTree entries: " << matchingHitsTree->GetEntries() << std::endl;
  matchingHitsTree->Print();

  if(matchingHitsTree->GetEntries() == 3 * nCubeHits) {
    isTreeFilled = true;
  }

}

TTree* FakeSFG::getMatchingHitsTree() {
  if(isTreeFilled) {
    return matchingHitsTree;
  } else {
    std::cerr<<"Error: Tree is empty. Returning a nullptr"<<std::endl;
    return nullptr;
  }
}