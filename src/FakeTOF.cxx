//
// Created by Vedantha Srinivas Kasturi on 14.04.2025.
//

#include "FakeTOF.hxx"

static const std::array<std::string, 6> PanelIdNames = {
        "TOP", "BOTTOM", "UPSTREAM", "DOWNSTREAM", "NORTH", "SOUTH"
};

std::string getPanelName(int id) {
    if (id >= 0 && id < static_cast<int>(PanelIdNames.size()))
        return PanelIdNames[id];
    else
        return "UNKNOWN";
}


FakeTOF::FakeTOF(std::string offset_type, double timeReso) : intrinsicTimeResolution(timeReso) {

    // Initialize the bars by fixing their centers
    std::vector<double> pos;
    // Set the first bar center positions
    pos = {TopBar0[0], TopBar0[1] , TopBar0[2]};
    Topbars[0] = pos;
    pos = {BottomBar0[0], BottomBar0[1] , BottomBar0[2]};
    Bottombars[0] = pos;
    pos = {UpstreamBar0[0], UpstreamBar0[1] , UpstreamBar0[2]};
    Upstreambars[0] = pos;
    pos = {DownstreamBar0[0], DownstreamBar0[1] , DownstreamBar0[2]};
    Downstreambars[0] = pos;
    pos = {NorthBar0[0], NorthBar0[1] , NorthBar0[2]};
    Northbars[0] = pos;
    pos = {SouthBar0[0], SouthBar0[1] , SouthBar0[2]};
    Southbars[0] = pos;

    // Set the rest of the bar center positions

    for (int i = 1; i < barsPerPanel; i++) {
        pos = {TopBar0[0], TopBar0[1] , TopBar0[2] + i * barSeparation};
        Topbars[i] = pos;

        pos = {BottomBar0[0], BottomBar0[1] , BottomBar0[2] + i * barSeparation};
        Bottombars[i] = pos;

        pos = {UpstreamBar0[0], UpstreamBar0[1] - i * barSeparation , UpstreamBar0[2]};
        Upstreambars[i] = pos;

        pos = {DownstreamBar0[0], DownstreamBar0[1] - i * barSeparation , DownstreamBar0[2]};
        Downstreambars[i] = pos;

        pos = {NorthBar0[0], NorthBar0[1] - i * barSeparation , NorthBar0[2]};
        Northbars[i] = pos;

        pos = {SouthBar0[0], SouthBar0[1] - i * barSeparation , SouthBar0[2]};
        Southbars[i] = pos;
    }

    std::cout << "Top bars initialized: " << Topbars.size() << std::endl;
    std::cout << "Bottom bars initialized: " << Bottombars.size() << std::endl;
    std::cout << "Upstream bars initialized: " << Upstreambars.size() << std::endl;
    std::cout << "Downstream bars initialized: " << Downstreambars.size() << std::endl;
    std::cout << "North bars initialized: " << Northbars.size() << std::endl;
    std::cout << "South bars initialized: " << Southbars.size() << std::endl;

    // For the TOF bars are the channels themselves

    // Initialize the tree
    std::cout << "Initializing tree...";
    initializeTree();
    std::cout << "Done!" << std::endl;

    rnd = TRandom3(42); // Seed the random number generator
    // Initialize time offsets
    setTimeOffsets(offset_type);

    std::cout <<"--------------------------------------\n Toy TOF initialized"<< std::endl;
    std::cout << " Bars per panel: " << barsPerPanel << std::endl;
    for(auto const& [key, val]: Topbars){
        std::cout << "Top bar " << key << " position: " << val[0] << ", " << val[1] << ", " << val[2] << std::endl;
    }
    for(auto const& [key, val]: Bottombars){
        std::cout << "Bottom bar " << key << " position: " << val[0] << ", " << val[1] << ", " << val[2] << std::endl;
    }
    for(auto const& [key, val]: Upstreambars){
        std::cout << "Upstream bar " << key << " position: " << val[0] << ", " << val[1] << ", " << val[2] << std::endl;
    }
    for(auto const& [key, val]: Downstreambars){
        std::cout << "Downstream bar " << key << " position: " << val[0] << ", " << val[1] << ", " << val[2] << std::endl;
    }
    for(auto const& [key, val]: Northbars){
        std::cout << "North bar " << key << " position: " << val[0] << ", " << val[1] << ", " << val[2] << std::endl;
    }
    for(auto const& [key, val]: Southbars){
        std::cout << "South bar " << key << " position: " << val[0] << ", " << val[1] << ", " << val[2] << std::endl;
    }
    std::cout << "Intrinsic time resolution: " << intrinsicTimeResolution << " ns" << std::endl;
    std::cout << "Time offsets initialized with " << offset_type << " method" <<std::endl;
    std::cout << "--------------------------------------\n";
    initializeConfigTree();
}

void FakeTOF::initializeTree() {
    matchingHitsTree = new TTree("matchingHitsTree", "Matching hits");
    matchingHitsTree->Branch("panel1", &thisMatchingHitsPair.panel1, "panel/I");
    matchingHitsTree->Branch("bar1", &thisMatchingHitsPair.bar1, "bar/I");
    matchingHitsTree->Branch("time1", &thisMatchingHitsPair.time1, "time/D");
    matchingHitsTree->Branch("panel2", &thisMatchingHitsPair.panel2, "panel/I");
    matchingHitsTree->Branch("bar2", &thisMatchingHitsPair.bar2, "bar/I");
    matchingHitsTree->Branch("time2", &thisMatchingHitsPair.time2, "time/D");
    matchingHitsTree->Branch("channelId1", &thisMatchingHitsPair.channel1, "channelId1/I");
    matchingHitsTree->Branch("channelId2", &thisMatchingHitsPair.channel2, "channelId2/I");
    matchingHitsTree->Branch("distance", &thisMatchingHitsPair.distance, "distance/D");
    matchingHitsTree->Branch("flightTimePostOffset", &thisMatchingHitsPair.flightTimePostOffset, "flightTimePostOffset/D");
    matchingHitsTree->Branch("flightTimeRaw", &thisMatchingHitsPair.flightTimeRaw, "flightTimeRaw/D");
    matchingHitsTree->Branch("correctedTime1", &thisMatchingHitsPair.correctedTime1, "correctedTime1/D");
    matchingHitsTree->Branch("correctedTime2", &thisMatchingHitsPair.correctedTime2, "correctedTime2/D");
    matchingHitsTree->Branch("hit1Position", &thisMatchingHitsPair.hit1Position, "hit1Position[3]/D");
    matchingHitsTree->Branch("hit2Position", &thisMatchingHitsPair.hit2Position, "hit2Position[3]/D");
    matchingHitsTree->Branch("DeltaTimePostOffset", &thisMatchingHitsPair.DeltaTimePostOffset, "DeltaTimePostOffset/D");
}

void FakeTOF::setTimeOffsets(std::string offset_type) {
    if(offset_type == "random") {
        setRandomTimeOffsets(1);
    } else if(offset_type == "gaussian") {
        setGaussianTimeOffsets(1);
    } else if(offset_type == "sawtooth") {
        setSawToothTimeOffsets(1);
    } else if(offset_type == "alternating") {
        setAlternatingTimeOffsets(1, 1);
    } else {
        std::cerr << "Error: Time offset type not recognized. Choose between random, gaussian, sawtooth, alternating" << std::endl;
    }

    // Push back the average to zero (we impose that the average of the offsets is 0)
    double average = 0;
    for(auto const& [key, val]: time_offsetsTop){
        average += val;
    }
    average /= time_offsetsTop.size();
    for(auto& [key, val]: time_offsetsTop){
        val -= average;
    }

    average = 0;
    for(auto const& [key, val]: time_offsetsBottom){
        average += val;
    }
    average /= time_offsetsBottom.size();
    for(auto& [key, val]: time_offsetsBottom){
        val -= average;
    }

    average = 0;
    for(auto const& [key, val]: time_offsetsUpstream){
        average += val;
    }

    average /= time_offsetsUpstream.size();
    for(auto& [key, val]: time_offsetsUpstream){
        val -= average;
    }

    average = 0;

    for(auto const& [key, val]: time_offsetsDownstream){
        average += val;
    }
    average /= time_offsetsDownstream.size();
    for(auto& [key, val]: time_offsetsDownstream){
        val -= average;
    }

    average = 0;

    for(auto const& [key, val]: time_offsetsNorth){
        average += val;
    }

    average /= time_offsetsNorth.size();

    for(auto& [key, val]: time_offsetsNorth){
        val -= average;
    }
    average = 0;

    for(auto const& [key, val]: time_offsetsSouth){
        average += val;
    }

    average /= time_offsetsSouth.size();
    for(auto& [key, val]: time_offsetsSouth){
        val -= average;
    }

    writableTimeOffsetsTop = new std::vector<double>;
    for(auto const& [key, val]: time_offsetsTop){
        writableTimeOffsetsTop->push_back(val);
    }

    writableTimeOffsetsBottom = new std::vector<double>;
    for(auto const& [key, val]: time_offsetsBottom){
        writableTimeOffsetsBottom->push_back(val);
    }

    writableTimeOffsetsUpstream = new std::vector<double>;
    for(auto const& [key, val]: time_offsetsUpstream){
        writableTimeOffsetsUpstream->push_back(val);
    }

    writableTimeOffsetsDownstream = new std::vector<double>;
    for(auto const& [key, val]: time_offsetsDownstream){
        writableTimeOffsetsDownstream->push_back(val);
    }

    writableTimeOffsetsNorth = new std::vector<double>;
    for(auto const& [key, val]: time_offsetsNorth){
        writableTimeOffsetsNorth->push_back(val);
    }

    writableTimeOffsetsSouth = new std::vector<double>;
    for(auto const& [key, val]: time_offsetsSouth){
        writableTimeOffsetsSouth->push_back(val);
    }

    writableTimeOffsets = new std::vector<double>;
    for(auto const& [key, val]: time_offsetsTop){
        writableTimeOffsets->push_back(val);
    }
    for(auto const& [key, val]: time_offsetsBottom){
        writableTimeOffsets->push_back(val);
    }
    for(auto const& [key, val]: time_offsetsUpstream){
        writableTimeOffsets->push_back(val);
    }

    for(auto const& [key, val]: time_offsetsDownstream){
        writableTimeOffsets->push_back(val);
    }
    for(auto const& [key, val]: time_offsetsNorth){
        writableTimeOffsets->push_back(val);
    }
    for(auto const& [key, val]: time_offsetsSouth){
        writableTimeOffsets->push_back(val);
    }

    std::cout<< "Time offsets Size" << writableTimeOffsets->size() << std::endl;

    std::cout << "Time offsets initialized" << std::endl;


}

void FakeTOF::setRandomTimeOffsets(double pitch) {

    for (auto const &[key, val]: Topbars) {
        if(time_offsetsTop.find(key) == time_offsetsTop.end() || time_offsetsTop[key] == 0 ) {
            time_offsetsTop[key] = rnd.Uniform(-pitch,pitch);
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Bottombars) {
        if(time_offsetsBottom.find(key) == time_offsetsBottom.end() || time_offsetsBottom[key] == 0 ) {
            time_offsetsBottom[key] = rnd.Uniform(-pitch,pitch);
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Upstreambars) {
        if(time_offsetsUpstream.find(key) == time_offsetsUpstream.end() || time_offsetsUpstream[key] == 0 ) {
            time_offsetsUpstream[key] = rnd.Uniform(-pitch,pitch);
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Downstreambars) {
        if(time_offsetsDownstream.find(key) == time_offsetsDownstream.end() || time_offsetsDownstream[key] == 0 ) {
            time_offsetsDownstream[key] = rnd.Uniform(-pitch,pitch);
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Northbars) {
        if(time_offsetsNorth.find(key) == time_offsetsNorth.end() || time_offsetsNorth[key] == 0 ) {
            time_offsetsNorth[key] = rnd.Uniform(-pitch,pitch);
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Southbars) {
        if (time_offsetsSouth.find(key) == time_offsetsSouth.end() || time_offsetsSouth[key] == 0) {
            time_offsetsSouth[key] = rnd.Uniform(-pitch, pitch);
        }
    }

}

void FakeTOF::setGaussianTimeOffsets(double sigma) {

    for (auto const &[key, val]: Topbars) {
        if(time_offsetsTop.find(key) == time_offsetsTop.end() || time_offsetsTop[key] == 0 ) {
            time_offsetsTop[key] = rnd.Gaus(0, sigma);
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Bottombars) {
        if(time_offsetsBottom.find(key) == time_offsetsBottom.end() || time_offsetsBottom[key] == 0 ) {
            time_offsetsBottom[key] = rnd.Gaus(0, sigma);
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Upstreambars) {
        if(time_offsetsUpstream.find(key) == time_offsetsUpstream.end() || time_offsetsUpstream[key] == 0 ) {
            time_offsetsUpstream[key] = rnd.Gaus(0, sigma);
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Downstreambars) {
        if(time_offsetsDownstream.find(key) == time_offsetsDownstream.end() || time_offsetsDownstream[key] == 0 ) {
            time_offsetsDownstream[key] = rnd.Gaus(0, sigma);
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Northbars) {
        if(time_offsetsNorth.find(key) == time_offsetsNorth.end() || time_offsetsNorth[key] == 0 ) {
            time_offsetsNorth[key] = rnd.Gaus(0, sigma);
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Southbars) {
        if (time_offsetsSouth.find(key) == time_offsetsSouth.end() || time_offsetsSouth[key] == 0) {
            time_offsetsSouth[key] = rnd.Gaus(0, sigma);
        } else {
            std::cerr << "Error: Time offset for channel " << key << " already set." << std::endl;
            return;
        }
    }

}

void FakeTOF::setSawToothTimeOffsets(double pitch) {
    int i= 0;
    for (auto const &[key, val]: Topbars) {
        if(time_offsetsTop.find(key) == time_offsetsTop.end() || time_offsetsTop[key] == 0 ) {
            time_offsetsTop[key] = (i%10)*pitch;
            i++;
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    i= 0;
    for (auto const &[key, val]: Bottombars) {
        if(time_offsetsBottom.find(key) == time_offsetsBottom.end() || time_offsetsBottom[key] == 0 ) {
            time_offsetsBottom[key] = (i%10)*pitch;
            i++;
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    i= 0;
    for (auto const &[key, val]: Upstreambars) {
        if(time_offsetsUpstream.find(key) == time_offsetsUpstream.end() || time_offsetsUpstream[key] == 0 ) {
            time_offsetsUpstream[key] = (i%10)*pitch;
            i++;
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    i= 0;
    for (auto const &[key, val]: Downstreambars) {
        if(time_offsetsDownstream.find(key) == time_offsetsDownstream.end() || time_offsetsDownstream[key] == 0 ) {
            time_offsetsDownstream[key] = (i%10)*pitch;
            i++;
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    i= 0;
    for (auto const &[key, val]: Northbars) {
        if(time_offsetsNorth.find(key) == time_offsetsNorth.end() || time_offsetsNorth[key] == 0 ) {
            time_offsetsNorth[key] = (i%10)*pitch;
            i++;
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    i= 0;
    for (auto const &[key, val]: Southbars) {
        if (time_offsetsSouth.find(key) == time_offsetsSouth.end() || time_offsetsSouth[key] == 0) {
            time_offsetsSouth[key] = (i%10)*pitch;
            i++;
        } else {
            std::cerr << "Error: Time offset for channel " << key << " already set." << std::endl;
            return;
        }
    }
}

void FakeTOF::setAlternatingTimeOffsets(double distance, double sigma) {
    for (auto const &[key, val]: Topbars) {
        if(time_offsetsTop.find(key) == time_offsetsTop.end() || time_offsetsTop[key] == 0 ) {
            double plusMinus = rnd.Uniform(-1,1);
            if(plusMinus > 0) {
                time_offsetsTop[key] = distance + rnd.Gaus(0, sigma);
            }else {
                time_offsetsTop[key] = -distance + rnd.Gaus(0, sigma);
            }
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Bottombars) {
        if(time_offsetsBottom.find(key) == time_offsetsBottom.end() || time_offsetsBottom[key] == 0 ) {
            double plusMinus = rnd.Uniform(-1,1);
            if(plusMinus > 0) {
                time_offsetsBottom[key] = distance + rnd.Gaus(0, sigma);
            }else {
                time_offsetsBottom[key] = -distance + rnd.Gaus(0, sigma);
            }
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Upstreambars) {
        if(time_offsetsUpstream.find(key) == time_offsetsUpstream.end() || time_offsetsUpstream[key] == 0 ) {
            double plusMinus = rnd.Uniform(-1,1);
            if(plusMinus > 0) {
                time_offsetsUpstream[key] = distance + rnd.Gaus(0, sigma);
            }else {
                time_offsetsUpstream[key] = -distance + rnd.Gaus(0, sigma);
            }
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Downstreambars) {
        if(time_offsetsDownstream.find(key) == time_offsetsDownstream.end() || time_offsetsDownstream[key] == 0 ) {
            double plusMinus = rnd.Uniform(-1,1);
            if(plusMinus > 0) {
                time_offsetsDownstream[key] = distance + rnd.Gaus(0, sigma);
            }else {
                time_offsetsDownstream[key] = -distance + rnd.Gaus(0, sigma);
            }
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Northbars) {
        if(time_offsetsNorth.find(key) == time_offsetsNorth.end() || time_offsetsNorth[key] == 0 ) {
            double plusMinus = rnd.Uniform(-1,1);
            if(plusMinus > 0) {
                time_offsetsNorth[key] = distance + rnd.Gaus(0, sigma);
            }else {
                time_offsetsNorth[key] = -distance + rnd.Gaus(0, sigma);
            }
        }else{
            std::cerr<<"Error: Time offset for channel "<<key<<" already set."<<std::endl;
            return;
        }
    }
    for (auto const &[key, val]: Southbars) {
        if (time_offsetsSouth.find(key) == time_offsetsSouth.end() || time_offsetsSouth[key] == 0) {
            double plusMinus = rnd.Uniform(-1, 1);
            if (plusMinus > 0) {
                time_offsetsSouth[key] = distance + rnd.Gaus(0, sigma);
            } else {
                time_offsetsSouth[key] = -distance + rnd.Gaus(0, sigma);
            }
        } else {
            std::cerr << "Error: Time offset for channel " << key << " already set." << std::endl;
            return;
        }
    }
}

FakeTOF::matchingHitsPair FakeTOF::simulateTimeHit() {
    if (verbose) {
        std::cout << "----------------------------------------\nDEBUG: simulating a time hit.\n";
    }
    // Get the random first bar
    int bar_f = rnd.Integer(barsPerPanel);
    int panel_f = rnd.Integer(6); // 0: Top, 1: Bottom, 2: Upstream, 3: Downstream, 4: North, 5: South

    // Get the random second bar
    int bar_s = rnd.Integer(barsPerPanel);
    int panel_s = rnd.Integer(6); // 0: Top, 1: Bottom, 2: Upstream, 3: Downstream, 4: North, 5: South

    // Get the time offset for the bar
    double time_offset_f = 0;
    double time_offset_s = 0;

    if (panel_f == PanelId::TOP) {
        time_offset_f = time_offsetsTop[bar_f];
    } else if (panel_f == PanelId::BOTTOM) {
        time_offset_f = time_offsetsBottom[bar_f];
    } else if (panel_f == PanelId::UPSTREAM) {
        time_offset_f = time_offsetsUpstream[bar_f];
    } else if (panel_f == PanelId::DOWNSTREAM) {
        time_offset_f = time_offsetsDownstream[bar_f];
    } else if (panel_f == PanelId::NORTH) {
        time_offset_f = time_offsetsNorth[bar_f];
    } else if (panel_f == PanelId::SOUTH) {
        time_offset_f = time_offsetsSouth[bar_f];
    }

    if (panel_s == PanelId::TOP) {
        time_offset_s = time_offsetsTop[bar_s];
    } else if (panel_s == PanelId::BOTTOM) {
        time_offset_s = time_offsetsBottom[bar_s];
    } else if (panel_s == PanelId::UPSTREAM) {
        time_offset_s = time_offsetsUpstream[bar_s];
    } else if (panel_s == PanelId::DOWNSTREAM) {
        time_offset_s = time_offsetsDownstream[bar_s];
    } else if (panel_s == PanelId::NORTH) {
        time_offset_s = time_offsetsNorth[bar_s];
    } else if (panel_s == PanelId::SOUTH) {
        time_offset_s = time_offsetsSouth[bar_s];
    }

    // Generate the first hit position
    // For upstream, downstream, top and bottom bars, the x coordinate can vary from -barlength/2 to barlength/2. Y and Z are fixed
    // For north and south bars, the Z coordinate can vary from -barlength/2 to barlength/2. X and Y are fixed

    double x_f = 0;
    double y_f = 0;
    double z_f = 0;

    if (panel_f == PanelId::TOP) {
        x_f = rnd.Uniform(-barLength / 2, barLength / 2);
        y_f = Topbars[bar_f][1];
        z_f = Topbars[bar_f][2];
    } else if( panel_f == PanelId::BOTTOM) {
        x_f = rnd.Uniform(-barLength / 2, barLength / 2);
        y_f = Bottombars[bar_f][1];
        z_f = Bottombars[bar_f][2];
    } else if (panel_f == PanelId::UPSTREAM) {
        x_f = rnd.Uniform(-barLength / 2, barLength / 2);
        y_f = Upstreambars[bar_f][1];
        z_f = Upstreambars[bar_f][2];
    } else if (panel_f == PanelId::DOWNSTREAM) {
        x_f = rnd.Uniform(-barLength / 2, barLength / 2);
        y_f = Downstreambars[bar_f][1];
        z_f = Downstreambars[bar_f][2];
    } else if (panel_f == PanelId::NORTH) {
        x_f = Northbars[bar_f][0];
        y_f = Northbars[bar_f][1];
        z_f = rnd.Uniform(-barLength / 2, barLength / 2);
    } else if (panel_f == PanelId::SOUTH) {
        x_f = Southbars[bar_f][0];
        y_f = Southbars[bar_f][1];
        z_f = rnd.Uniform(-barLength / 2, barLength / 2);
    }

    // Generate the second hit position

    double x_s = 0;
    double y_s = 0;
    double z_s = 0;

    if (panel_s == PanelId::TOP) {
        x_s = rnd.Uniform(-barLength / 2, barLength / 2);
        y_s = Topbars[bar_s][1];
        z_s = Topbars[bar_s][2];
    } else if( panel_s == PanelId::BOTTOM) {
        x_s = rnd.Uniform(-barLength / 2, barLength / 2);
        y_s = Bottombars[bar_s][1];
        z_s = Bottombars[bar_s][2];
    } else if (panel_s == PanelId::UPSTREAM) {
        x_s = rnd.Uniform(-barLength / 2, barLength / 2);
        y_s = Upstreambars[bar_s][1];
        z_s = Upstreambars[bar_s][2];
    } else if (panel_s == PanelId::DOWNSTREAM) {
        x_s = rnd.Uniform(-barLength / 2, barLength / 2);
        y_s = Downstreambars[bar_s][1];
        z_s = Downstreambars[bar_s][2];
    } else if (panel_s == PanelId::NORTH) {
        x_s = Northbars[bar_s][0];
        y_s = Northbars[bar_s][1];
        z_s = rnd.Uniform(-barLength / 2, barLength / 2);
    } else if (panel_s == PanelId::SOUTH) {
        x_s = Southbars[bar_s][0];
        y_s = Southbars[bar_s][1];
        z_s = rnd.Uniform(-barLength / 2, barLength / 2);
    }

    if(panel_s == panel_f){
        return simulateTimeHit();
    }

    // Calculate the distance between the two hits

    double distance = sqrt(pow(x_f - x_s, 2) + pow(y_f - y_s, 2) + pow(z_f - z_s, 2));

    if(verbose) {
        if(distance < 100){
            if(panel_f == PanelId::TOP && panel_s == PanelId::BOTTOM) std::cout << "DEBUG: Top and Bottom bars hit. Distance too small: " << distance << std::endl;
            if(panel_f == PanelId::UPSTREAM && panel_s == PanelId::DOWNSTREAM) std::cout << "DEBUG: Upstream and Downstream bars hit. Distance too small: " << distance << std::endl;
            if(panel_f == PanelId::NORTH && panel_s == PanelId::SOUTH) std::cout << "DEBUG: North and South bars hit. Distance too small: " << distance << std::endl;
        }

        if(distance > (sqrt(3)*barLength)) {
            std::cout << "DEBUG: Distance too large: " << distance << std::endl;
        }
    }

    // Calculate the time of flight
    double time_Raw = distance / speed_of_light;

    matchingHitsPair pair;

    pair.panel1 = panel_f;
    pair.bar1 = bar_f;
    pair.time1 = 0 + time_offset_f + rnd.Gaus(0, intrinsicTimeResolution);
    pair.panel2 = panel_s;
    pair.bar2 = bar_s;
    pair.time2 = time_Raw + time_offset_s + rnd.Gaus(0, intrinsicTimeResolution);
    pair.channel1 = getChannelID(panel_f, bar_f); // Get the channel ID for the first hit
    pair.channel2 = getChannelID(panel_s, bar_s); // Get the channel ID for the second hit
    pair.distance = distance;
    pair.flightTimePostOffset = pair.time2 - pair.time1;
    pair.flightTimeRaw = time_Raw;
    pair.correctedTime1 = pair.time1 - time_offset_f;
    pair.correctedTime2 = pair.time2 - time_offset_s;
    pair.hit1Position[0] = x_f;
    pair.hit1Position[1] = y_f;
    pair.hit1Position[2] = z_f;
    pair.hit2Position[0] = x_s;
    pair.hit2Position[1] = y_s;
    pair.hit2Position[2] = z_s;
    pair.DeltaTimePostOffset = pair.time2 - pair.time1 - distance / speed_of_light;

    if(verbose) {
        std::cout << "DEBUG: panel1: " << getPanelName(pair.panel1) << ", bar1: " << pair.bar1 << ", Channel Id: "<< pair.channel1 << ", time1: " << pair.time1 << std::endl;
        std::cout << "DEBUG: panel2: " << getPanelName(pair.panel2) << ", bar2: " << pair.bar2 << ", Channel Id: "<< pair.channel2 << ", time2: " << pair.time2 << std::endl;
        std::cout << "DEBUG: distance: " << pair.distance << std::endl;
        std::cout << "DEBUG: hit1Position: " << pair.hit1Position[0] << ", " << pair.hit1Position[1] << ", " << pair.hit1Position[2] << std::endl;
        std::cout << "DEBUG: hit2Position: " << pair.hit2Position[0] << ", " << pair.hit2Position[1] << ", " << pair.hit2Position[2] << std::endl;
    }

    return pair;
}


void FakeTOF::generateMatchingHitsTree(int NHits) {

    if (!matchingHitsTree){
        std::cerr << "Error: Tree not initialized. Please initialize the tree before generating matching hits." << std::endl;
        return;
    }

    std::cout << "Generating matching hits pairs for " << NHits << " hits..." << std::endl;
    int progress = 0;
    for (int i = 0; i < NHits; i++) {
        if (i % (NHits / 100) == 0) {
            progress++;
            std::string bar = std::string(progress, '|');
            std::cout << "\r" << bar << " " <<progress<<" %" <<std::flush;
        }

        thisMatchingHitsPair = simulateTimeHit();
        if (verbose) {
            std::cout << "Panel1: " << getPanelName(thisMatchingHitsPair.panel1) << ", Bar1: " << thisMatchingHitsPair.bar1 << ", ChannelID: "<< thisMatchingHitsPair.channel1 << ", Time1: " << thisMatchingHitsPair.time1 << std::endl;
            std::cout << "Panel2: " << getPanelName(thisMatchingHitsPair.panel2) << ", Bar2: " << thisMatchingHitsPair.bar2 << ", ChannelID: "<< thisMatchingHitsPair.channel2 << ", Time2: " << thisMatchingHitsPair.time2 << std::endl;
            std::cout << "Distance: " << thisMatchingHitsPair.distance << std::endl;
        }
        matchingHitsTree->Fill();
    }
    std::cout << "\rDONE!" << std::endl;
    std::cout << "\n Simulated " << NHits << " hits." << std::endl;
    std::cout << "Tree filled with " << matchingHitsTree->GetEntries() << " entries." << std::endl;
    matchingHitsTree->Print(); // Print the tree structure
}

void FakeTOF::initializeConfigTree() {

    configTree = new TTree("configTree", "Tree containing the configuration of the bars in the Toy TOF");

    configTree->Branch("Topbars", &Topbars, "Topbars[10][3]/D");
    configTree->Branch("Bottombars", &Bottombars, "Bottombars[10][3]/D");
    configTree->Branch("Upstreambars", &Upstreambars, "Upstreambars[10][3]/D");
    configTree->Branch("Downstreambars", &Downstreambars, "Downstreambars[10][3]/D");
    configTree->Branch("Northbars", &Northbars, "Northbars[10][3]/D");
    configTree->Branch("Southbars", &Southbars, "Southbars[10][3]/D");
    configTree->Branch("intrinsicTimeResolution", &intrinsicTimeResolution, "intrinsicTimeResolution/D");
    configTree->Branch("speed_of_light", &speed_of_light, "speed_of_light/D");
    configTree->Branch("barLength", &barLength, "barLength/D");
    configTree->Branch("barSeparation", &barSeparation, "barSeparation/D");
    configTree->Branch("barsPerPanel", &barsPerPanel, "barsPerPanel/I");
    configTree->Branch("TopBar0", &TopBar0, "TopBar0[3]/D");
    configTree->Branch("BottomBar0", &BottomBar0, "BottomBar0[3]/D");
    configTree->Branch("UpstreamBar0", &UpstreamBar0, "UpstreamBar0[3]/D");
    configTree->Branch("DownstreamBar0", &DownstreamBar0, "DownstreamBar0[3]/D");
    configTree->Branch("NorthBar0", &NorthBar0, "NorthBar0[3]/D");
    configTree->Branch("SouthBar0", &SouthBar0, "SouthBar0[3]/D");
    configTree->Branch("time_offsetsTop", &writableTimeOffsetsTop);
    configTree->Branch("time_offsetsBottom", &writableTimeOffsetsBottom);
    configTree->Branch("time_offsetsUpstream", &writableTimeOffsetsUpstream);
    configTree->Branch("time_offsetsDownstream", &writableTimeOffsetsDownstream);
    configTree->Branch("time_offsetsNorth", &writableTimeOffsetsNorth);
    configTree->Branch("time_offsetsSouth", &writableTimeOffsetsSouth);
    configTree->Branch("time_offsets", &writableTimeOffsets);

    configTree->Fill();
}

