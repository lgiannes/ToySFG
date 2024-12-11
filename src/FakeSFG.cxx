//
// Created by Lorenzo Giannessi on 26.11.2024.
//

#include "FakeSFG.hxx"


FakeSFG::FakeSFG(int n_x, int n_y, int n_z, std::string offset_type, double timeReso, double timeStep) : n_x(n_x), n_y(n_y), n_z(n_z), intrinsicTimeResolution(timeReso), timeStep(timeStep){

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
    nCubes = cubes.size();
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
    nChannels = channels.size();

    std::cout<<"Initializing tree...";
    initializeTree();
    isTreeFilled = false;
    std::cout<<" done"<<std::endl;

    if(timeStep < 0){
      isTimeDiscrete = false;
    }else{
      isTimeDiscrete = true;
    }
    // Initialize random generator
    rnd = TRandom3(42);
    // Initialize time offsets
    setTimeOffsets(offset_type);
    // Generate map cube to channels
    generateMapCubeToChannels();

    std::cout<<"-------------------------------------\nToy SFGD initialized"<<std::endl;
    std::cout<<"n_x: "<<n_x<<" n_y: "<<n_y<<" n_z: "<<n_z<<std::endl;
    std::cout<<"Time resolution: "<<intrinsicTimeResolution<<" ns"<<std::endl;
    if(isTimeDiscrete){
      std::cout<<"Time is discretized with time step: "<<timeStep<<" ns"<<std::endl;
    }else{
      std::cout<<"Time is not discretized."<<std::endl;
    }
    std::cout<<"Time offsets initialized with "<<offset_type<<" method"<<std::endl;
    std::cout<<"-------------------------------------\n";
    // Initialize config tree
    initializeConfigTree();


}

void FakeSFG::initializeConfigTree() {
  configTree = new TTree("configTree", "Tree containing the configuration of the Toy SFGD");
  configTree->Branch("n_x", &n_x, "n_x/I");
  configTree->Branch("n_y", &n_y, "n_y/I");
  configTree->Branch("n_z", &n_z, "n_z/I");
  configTree->Branch("nChannels", &nChannels, "nChannels/I");
  configTree->Branch("nCubes", &nCubes, "nCubes/I");
  configTree->Branch("timeResolution", &intrinsicTimeResolution, "timeResolution/D");
  configTree->Branch("timeStep", &timeStep, "timeStep/D");
  if (!writableTimeOffsets){
    std::cerr<<"Error: Time offsets not initialized"<<std::endl;
    throw std::runtime_error("Time offsets not initialized! Cannot save config tree");
  }
  configTree->Branch("timeOffsets", &writableTimeOffsets);
  configTree->Fill();
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

void FakeSFG::resetTree() {
    matchingHitsTree->Reset();
    isTreeFilled = false;
}

void FakeSFG::setTimeOffsets(std::string offset_type) {
  if(offset_type == "random") {
    setRandomTimeOffsets(1);
  } else if(offset_type == "gaussian") {
    setGaussianTimeOffsets(1);
  } else if(offset_type == "sawtooth") {
    setSawToothTimeOffsets(1);
  } else if(offset_type == "alternating") {
    setAlternatingTimeOffsets(1, 1);
  } else {
    std::cerr<<"Error: Time offset type not recognized. Choose between random, gaussian, sawtooth, alternating"<<std::endl;
  }
  writableTimeOffsets = new std::vector<double>;
  for(auto const& [key, val]: time_offsets){
    writableTimeOffsets->push_back(val);
  }
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

    channels_reading = mapCubeToChannels[random];

    std::pair< unsigned int , std::vector<unsigned int> > result = std::make_pair(random, channels_reading); // cube number, channels reading the cube


    return result;

}


void FakeSFG::generateMapCubeToChannels() {
  std::cout<<"Generating map cube to channels..."<<std::endl;
  // Progress bar
  int progress = 0;
  for (auto const& [key, val] : cubes) {
      if(progress % (nCubes/100) == 0) {
        std::string bar = std::string(progress/(nCubes/100), '|');
        std::cout << "\r" << bar << " " <<progress/(nCubes/100)<<" %" <<std::flush;
      }
      progress++;
      std::vector<unsigned int> channels_reading;
      for (auto const& [key2, val2] : channels) {
          std::vector<int> pos_channel = val2;
          if (channels_projection[key2] == X && pos_channel[1] == val[1] && pos_channel[2] == val[2]) {
              channels_reading.push_back(key2);
          }
          if (channels_projection[key2] == Y && pos_channel[2] == val[2] && pos_channel[0] == val[0]) {
              channels_reading.push_back(key2);
          }
          if (channels_projection[key2] == Z && pos_channel[1] == val[1] && pos_channel[0] == val[0]) {
              channels_reading.push_back(key2);
          }
      }
      if(channels_reading.size() != 3){
          std::cerr<<"ERROR: Cube "<<key<<" is not read by 3 channels"<<std::endl;
      }
      mapCubeToChannels[key] = channels_reading;
  }
  std::cout<<"Map cube to channels generated"<<std::endl;
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
  times[0] = time_x + rnd.Gaus(0, intrinsicTimeResolution) + time_offsets[channels_reading[0]];
  double time_y = getDistance(cube, channels_reading[1]) * cubeWidth / speed_of_light;
  times[1] = time_y + rnd.Gaus(0, intrinsicTimeResolution) + time_offsets[channels_reading[1]];
  double time_z = getDistance(cube, channels_reading[2]) * cubeWidth / speed_of_light;
  times[2] = time_z + rnd.Gaus(0, intrinsicTimeResolution) + time_offsets[channels_reading[2]];

  if(isTimeDiscrete){
    for(int i=0;i<3;i++){
      times[i] = round(times[i]/timeStep)*timeStep;
    }
  }

  std::vector<matchingHitsPair> ThreeHitPairs;
  for(int i = 0; i < 3; i++) {
    matchingHitsPair pair;
    pair.channel1 = channels_reading[i];
    pair.time1 = times[i];
    pair.distance1 = getDistance(cube, channels_reading[i])*cubeWidth;
    pair.channel2 = channels_reading[(i + 1) % 3];
    pair.time2 = times[(i + 1) % 3];
    pair.distance2 = getDistance(cube, channels_reading[(i + 1) % 3])*cubeWidth;
    pair.projection1 = channels_projection[channels_reading[i]];
    pair.projection2 = channels_projection[channels_reading[(i + 1) % 3]];
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

  std::cout<<" Generating matching hits pairs for "<<nCubeHits<<" cube hits..."<<std::endl;
  int progress = 0;
  for(int i = 0; i < nCubeHits; i++) {
    // Progress bar that prints a | every % of the total hits
      if(i % (nCubeHits/100) == 0) {
        progress++;
        std::string bar = std::string(progress, '|');
        std::cout << "\r" << bar << " " <<progress<<" %" <<std::flush;
      }
    std::vector<matchingHitsPair> hits = simulateTimeHit();
    for (auto const &hit: hits) {
      thisMatchingHitsPair = hit;
      if(verbose>0){
        std::cout<<"time1: "<<thisMatchingHitsPair.time1<<" "<<"time2: "<<thisMatchingHitsPair.time2<<std::endl;
        std::cout<<"distance1: "<<thisMatchingHitsPair.distance1<<" "<<"distance2: "<<thisMatchingHitsPair.distance2<<std::endl;
        std::cout<<"channel1: "<<thisMatchingHitsPair.channel1<<" "<<"channel2: "<<thisMatchingHitsPair.channel2<<std::endl;
        std::cout<<"projection1: "<<thisMatchingHitsPair.projection1<<" "<<"projection2: "<<thisMatchingHitsPair.projection2<<std::endl;
      }
      matchingHitsTree->Fill();
    }
  }
  std::cout << "\rDONE!" << std::endl;

  std::cout << "\nSimulated " << nCubeHits << " cube hits.\n";
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

void FakeSFG::setIntrinsicTimeResolution(double res) {
  intrinsicTimeResolution = res;
}

void FakeSFG::setTimeStep(double ts) {
  if(ts<=0){
    timeStep = 0;
    isTimeDiscrete = false;
  }else{
    isTimeDiscrete = true;
    timeStep = ts;
  }
}

void FakeSFG::setRandomTimeOffsets(double pitch) {
  for (auto const &[key, val]: channels) {
    if(time_offsets.find(key) == time_offsets.end() || time_offsets[key] == 0 ) {
      time_offsets[key] = rnd.Uniform(-pitch,pitch);
    }else{
      std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
      return;
    }
  }
}

void FakeSFG::setGaussianTimeOffsets(double sigma) {
  for (auto const &[key, val]: channels) {
    if(time_offsets.find(key) == time_offsets.end() || time_offsets[key] == 0 ) {
      time_offsets[key] = rnd.Gaus(0, sigma);
    }else{
      std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
      return;
    }
  }
}

void FakeSFG::setSawToothTimeOffsets(double pitch) {
  int i= 0;
  for (auto const &[key, val]: channels) {
    if(time_offsets.find(key) == time_offsets.end() || time_offsets[key] == 0 ) {
      time_offsets[key] = (i%10)*pitch;
      i++;
    }else{
      std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
      return;
    }
  }
}

void FakeSFG::setAlternatingTimeOffsets(double distance, double sigma) {
  for (auto const &[key, val]: channels) {
    if(time_offsets.find(key) == time_offsets.end() || time_offsets[key] == 0 ) {
      double plusMinus = rnd.Uniform(-1,1);
      if(plusMinus > 0) {
        time_offsets[key] = distance + rnd.Gaus(0, sigma);
      }else {
        time_offsets[key] = -distance + rnd.Gaus(0, sigma);
      }
    }else{
      std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
      return;
    }
  }
}

void FakeSFG::resetTimeOffsets() {
  for (auto const &[key, val]: channels) {
    time_offsets[key] = 0;
  }
}