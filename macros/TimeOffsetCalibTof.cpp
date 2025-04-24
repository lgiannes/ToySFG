//
// Created by Vedantha Srinivas Kasturi on 22.04.2025.
//

void getPanelAndBar(int channelID, int barsPerPanel, int &panel, int &bar) {
if (channelID < 0 || channelID >= 6 * barsPerPanel) {
panel = -1;
bar = -1;
return;
}
panel = channelID / barsPerPanel;
bar = channelID % barsPerPanel;
}


double TimeOffsetCalibTof(std::string filename = "", bool saveImage = true, bool batchMode = false, int maxEntries = -1) {
    if(batchMode){
        gROOT->SetBatch(kTRUE);
    }

    TFile *file = new TFile(filename.c_str(), "READ");
    TTree *tree = (TTree *) file->Get("matchingHitsTree");  

    int panel1;
    int panel2;
    int bar1;
    int bar2;
    double time1;
    double time2;
    int channel1;
    int channel2;
    double distance;
    double correctedTime1;
    double correctedTime2;
    double DeltaTimePreComp;

    tree->SetBranchAddress("panel1", &panel1);
    tree->SetBranchAddress("panel2", &panel2);
    tree->SetBranchAddress("bar1", &bar1);
    tree->SetBranchAddress("bar2", &bar2);
    tree->SetBranchAddress("time1", &time1);
    tree->SetBranchAddress("time2", &time2);
    tree->SetBranchAddress("channelId1", &channel1);
    tree->SetBranchAddress("channelId2", &channel2);
    tree->SetBranchAddress("distance", &distance);
    tree->SetBranchAddress("correctedTime1", &correctedTime1);
    tree->SetBranchAddress("correctedTime2", &correctedTime2);
    tree->SetBranchAddress("DeltaTimePostOffset", &DeltaTimePreComp);

    int nEntries = tree->GetEntries();

    if(maxEntries > 0 && nEntries > maxEntries){
        nEntries = maxEntries;
    }

    double v = 30; // cm/ns

    // open config tree

    TTree *configTree = (TTree *) file->Get("configTree");
    vector<double> *writableTimeOffsets = 0;
    double timeResolution,speedOfLight;
    int barsPerPanel;

    configTree->SetBranchAddress("intrinsicTimeResolution", &timeResolution);
    configTree->SetBranchAddress("barsPerPanel", &barsPerPanel);
    configTree->SetBranchAddress("time_offsets", &writableTimeOffsets);
    configTree->SetBranchAddress("speed_of_light", &speedOfLight);
    configTree->GetEntry(0);

    map<int, double> inputTimeOffsets;
    for (int i = 0; i < writableTimeOffsets->size(); i++) {
        inputTimeOffsets[i] = writableTimeOffsets->at(i);
        cout << "Channel " << i << ": " << inputTimeOffsets[i] << endl;
    }
    int nChannels = inputTimeOffsets.size();

    // Start time offset calibration
    std::map<unsigned int, double> ch_sumDelta; // sum of delta
    std::map<unsigned int, double> ch_squaredDelta; // squared
    std::map<unsigned int, int> ch_entries; // entries
    std::map<unsigned int, int> ch_entries_persistent; // entries. to save after the iteration loop
    std::map<unsigned int, double> ch_offset; // sum of delta/N
    std::map<unsigned int, double> ch_offsetError; // squared/N
    std::map<unsigned int, double> lookUpTable; // final lookup table
    std::map<unsigned int, double> single_ch_offsetError; // error of each channel
    std::map<unsigned int, double> ch_sumDelta_t1; // sum of delta
    std::map<unsigned int, double> ch_sumDelta_t2; // sum of delta
    std::map<unsigned int, double> ch_sumDelta_t1_entries; // entries
    std::map<unsigned int, double> ch_sumDelta_t2_entries; // entries
    std::map<unsigned int, double> ch_sumDelta_ratio; // ratio of delta

    int maxIteration = 100;
    int iter = 0;
    double alpha = 0.05;
    int usedChannels = 0;

    int testchannel = 10;
    TGraph* g_computed_offset_testchannel = new TGraph();
    TGraph* g_delta_t1 = new TGraph();
    TGraph* g_delta_t2 = new TGraph();
    TGraph* g_delta_t1t2ratio = new TGraph();
    TGraph* g_magnitudeDeltaTk = new TGraph();
    TGraph* g_magnitudeDeltaTk_ratio = new TGraph();
    TH1F* h_offset_correction_ratio_last = new TH1F("h_offset_correction_ratio_last", "Offset correction ratio last iteration", 100, 0., 1);
    vector<double> offset_correction_ratio_secondlast;

    cout << "Start time offset calibration" << endl;
    cout << "Matching hits in dataset: " << nEntries << endl;
    cout << "Toy TOF channels: " << nChannels << endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // START OF ITERATION LOOP
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    while (iter < maxIteration){
        for (int hit_i = 0; hit_i < nEntries; hit_i++) {
            tree->GetEntry(hit_i);

            time1 -= lookUpTable[channel1];
            time2 -= lookUpTable[channel2];

            double delta_t_1 = time1 - correctedTime1;
            double delta_t_2 = time2 - correctedTime2;
            delta_t_1 *= alpha;
            delta_t_2 *= alpha;

            double delta = time2 - time1- (distance) / speedOfLight;
            if(iter == 0) {
                if (delta != DeltaTimePreComp) {
                    cout << "Delta is not equal to DeltaTimePreComp" << endl;
                    cout << "Delta: " << delta << " DeltaTime Pre Computed: " << DeltaTimePreComp << endl;
                }
            }
            delta *= alpha;

            ch_sumDelta_t1[channel1] += delta_t_1;
            ch_sumDelta_t2[channel1] += delta_t_2;
            ch_sumDelta_ratio[channel1] += delta_t_2 / delta_t_1;
            ch_sumDelta_t1_entries[channel1]++;
            ch_sumDelta_t2_entries[channel1]++;

            ch_sumDelta[channel1] -= delta;
            ch_squaredDelta[channel1] += delta * delta;
            ch_entries[channel1]++;
            ch_sumDelta_t1[channel2] += delta_t_2;
            ch_sumDelta_t2[channel2] += delta_t_1;
            ch_sumDelta_t1_entries[channel2]++;
            ch_sumDelta_t2_entries[channel2]++;
            ch_sumDelta[channel2] += delta;
            ch_squaredDelta[channel2] += delta * delta;
            ch_entries[channel2]++;
        } // end of loop over hits

        usedChannels = ch_entries.size();

        for (auto& entry : ch_sumDelta) {
            entry.second /= ch_entries[entry.first];
            ch_offset[entry.first] = entry.second;
            // Standard deviation is sqrt(E[X^2] - E[X]^2)
            ch_offsetError[entry.first] = sqrt(ch_squaredDelta[entry.first] / ch_entries[entry.first] - entry.second * entry.second);
            //std::cout << "nEntries is: " << ch_entries[entry.first] << std::endl;
            ch_sumDelta_t1[entry.first] /= ch_sumDelta_t1_entries[entry.first];
            ch_sumDelta_t2[entry.first] /= ch_sumDelta_t2_entries[entry.first];
            ch_sumDelta_ratio[entry.first] /= ch_sumDelta_t1_entries[entry.first];
        }
        // Compute magnitude of DeltaT^k vector
        double magDeltaTk = 0;
        for(auto deltaTk : ch_sumDelta_t1){
            magDeltaTk += deltaTk.second*deltaTk.second;
        }
        magDeltaTk = sqrt(magDeltaTk);
        cout << "|Delta T^{k}| = " << magDeltaTk << endl;
        cout << "average Delta1: " << ch_sumDelta_t1[testchannel] << " average Delta2: " << ch_sumDelta_t2[testchannel] << endl;
        cout << "average Delta: " << ch_offset[testchannel] << endl;
        cout << "Total correction: " << lookUpTable[testchannel] << " actual offset: " << inputTimeOffsets[testchannel] << endl;
        cout << "End of iteration " << iter << ", used " << usedChannels << " channels." << endl;

        g_magnitudeDeltaTk->SetPoint(iter,iter, magDeltaTk);
        if(iter>0) g_magnitudeDeltaTk_ratio->SetPoint(iter-1,iter-1,magDeltaTk / g_magnitudeDeltaTk->GetY()[iter-1]);
        g_delta_t1->SetPoint(iter, iter, ch_sumDelta_t1[testchannel]);
        g_delta_t2->SetPoint(iter, iter, ch_sumDelta_t2[testchannel]);
        g_delta_t1t2ratio->SetPoint(iter, iter, ch_sumDelta_ratio[testchannel]);

        // Shift all the offsets to have average = 0
        double sum = 0;
        for (auto& entry : ch_offset) {
            sum += entry.second;
        }
        double average = sum/usedChannels;
        for (auto& entry : ch_offset) {
            entry.second -= average;
        }

        // Fill lookup table
        for (auto& entry : ch_offset) {
//      std::cout << "Entries of channel " << entry.first << " is: " << ch_entries[entry.first] << std::endl;
            lookUpTable[entry.first] += entry.second;
            if (iter == maxIteration-2){
                offset_correction_ratio_secondlast.push_back(entry.second);
            }
            if (iter == maxIteration-1){
                h_offset_correction_ratio_last->Fill(entry.second/offset_correction_ratio_secondlast[entry.first]);
            }
        }

        g_computed_offset_testchannel->SetPoint(iter, iter, lookUpTable[testchannel]);

        if (iter != maxIteration-1){
            ch_sumDelta.clear();
            ch_squaredDelta.clear();
            ch_entries.clear();
            ch_sumDelta_t1.clear();
            ch_sumDelta_t2.clear();
            ch_sumDelta_t1_entries.clear();
            ch_sumDelta_t2_entries.clear();
        }
        iter++;
    } // end of iteration loop

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // END OF ITERATION LOOP
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if (usedChannels < nChannels){
        int notCalibrated = nChannels - usedChannels;
        cout << "Not calibrated channels: " << notCalibrated << endl;
        for(int i_ch = 0; i_ch < nChannels; i_ch++){
            if (ch_entries.find(i_ch) == ch_entries.end()){
                cout << i_ch << " ";
            }
        }
        cout << endl;
    }

    TH1F* h_offset_difference = new TH1F("h_offset_difference", "Difference between input and output offsets; #Delta T^{0}_{set} - #Delta T^{0}_{comp}  [ns]", 5000, -1, 1);
    TH1F* h_inputOffsets = new TH1F("h_inputOffsets", "Input offsets; #Delta T^{0} [ns]", 25, -3, 3);
    TH1F* h_outputOffsets = new TH1F("h_outputOffsets", "Output offsets; #Delta T^{0} [ns]", 25, -3, 3);
    TGraph* g_inputOffsets = new TGraph();
    TGraph* g_outputOffsets = new TGraph();
    TGraph* g_residuals = new TGraph();

    // Print out lookup table and compare to the input
    for(int channel = 0; channel < nChannels; channel++){
//    cout<<"Channel "<<channel<<": "<<lookUpTable[channel]<<" input: "<<inputTimeOffsets[channel]<<""<<endl;
        h_offset_difference->Fill(lookUpTable[channel] - inputTimeOffsets[channel]);
        h_inputOffsets->Fill(inputTimeOffsets[channel]);
        h_outputOffsets->Fill(lookUpTable[channel]);
        g_inputOffsets->SetPoint(g_inputOffsets->GetN(), channel, inputTimeOffsets[channel]);
        g_outputOffsets->SetPoint(g_outputOffsets->GetN(), channel, lookUpTable[channel]);
        g_residuals->SetPoint(g_residuals->GetN(), ch_entries[channel], 1000*(lookUpTable[channel] - inputTimeOffsets[channel]) );
    }

    TH1D* h_deltaBeforeCalib = new TH1D("h_deltaBeforeCalib", "Delta before calibration", 420, -10, 10);
    TH1D* h_deltaAfterCalib = new TH1D("h_deltaAfterCalib", "Delta after calibration", 420, -10, 10);
    for (int hit_i = 0; hit_i<nEntries; hit_i++) {
        tree->GetEntry(hit_i);
        double delta = time2 - time1 - (distance) / speedOfLight;
        double correctedDelta = time2 - time1 - (distance) / speedOfLight - (lookUpTable[channel2] - lookUpTable[channel1]);
        h_deltaBeforeCalib->Fill(delta);
        h_deltaAfterCalib->Fill(correctedDelta);
    }

    TCanvas* c_offDiff = new TCanvas("c_offDiff", "c_offDiff", 1800, 600);
    c_offDiff->Divide(2,1);
    c_offDiff->cd(1);
    h_offset_difference->Draw();
    h_offset_difference->Fit("gaus","Q");
//  double sigma_offset = h_offset_difference->GetFunction("gaus")->GetParameter(2);
    double sigma_offset = h_offset_difference->GetRMS();
    cout << "Mean of Offset_{measured} - Offset_{set} = " <<1000*h_offset_difference->GetMean() << " ps." << endl;
    ofstream ofs;
//  ofs.open("Reso_hits_sigmaOffsets.txt", std::ofstream::out | std::ofstream::app);
//  ofs << timeResolution << " " << nEntries << " " << sigma_offset << endl;
    c_offDiff->cd(2);
    h_deltaAfterCalib->SetLineColor(kBlue);
    h_deltaAfterCalib->Draw();
    h_deltaAfterCalib->Fit("gaus","Q");
    h_deltaBeforeCalib->SetLineColor(kRed);
    h_deltaBeforeCalib->Draw("same");
    h_deltaBeforeCalib->Fit("gaus","Q");
    double sigma_before = h_deltaBeforeCalib->GetFunction("gaus")->GetParameter(2);
    double sigma_after = h_deltaAfterCalib->GetFunction("gaus")->GetParameter(2);
    TLatex latex;
    latex.SetTextSize(0.03);
    latex.SetNDC();
    latex.DrawLatex(0.2, 0.91, Form("Before calib sigma: %.2f", sigma_before));
    latex.DrawLatex(0.2, 0.85, Form("After calib sigma: %.2f", sigma_after));
    TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
    leg->AddEntry(h_deltaBeforeCalib, "Before calibration", "l");
    leg->AddEntry(h_deltaAfterCalib, "After calibration", "l");
    leg->Draw();
    string reducedFilename = filename.substr(filename.find_last_of("/")+1);
    h_deltaAfterCalib->SetTitle(reducedFilename.c_str());
    h_deltaBeforeCalib->SetTitle(reducedFilename.c_str());

    if (saveImage)
        c_offDiff->SaveAs(Form("TimeOffsetCalib_%s.png", reducedFilename.c_str()));

    TCanvas* c_offset = new TCanvas("c_offset", "c_offset", 1600, 600);
    c_offset->Divide(2,1);
    c_offset->cd(1);
    g_computed_offset_testchannel->SetMarkerStyle(20);
    g_computed_offset_testchannel->SetMarkerSize(0.5);
    g_computed_offset_testchannel->SetMarkerColor(1);
    g_computed_offset_testchannel->SetLineColor(1);
    g_computed_offset_testchannel->SetTitle(Form("Offset of channel %d;Iteration;Offset [ns]", testchannel));
    g_computed_offset_testchannel->Draw("APL");
    TLine* line = new TLine(0, inputTimeOffsets[testchannel], maxIteration+1, inputTimeOffsets[testchannel]);
    line->SetLineColor(kRed);
    line->Draw();
    // Print difference between input and output (last iteration)
    double diff = lookUpTable[testchannel] - inputTimeOffsets[testchannel];
    TLatex latex2;
    latex2.SetTextSize(0.03);
    latex2.SetNDC();
    latex2.DrawLatex(0.2, 0.91, Form("Difference between input and output: %.1f ps", 1000*diff));
    c_offset->cd(2);
    TGraph* g_offset_correction = new TGraph();
    for (int i = 1; i < g_computed_offset_testchannel->GetN(); i++) {
        g_offset_correction->SetPoint(g_offset_correction->GetN(), i, abs(g_computed_offset_testchannel->GetY()[i] - g_computed_offset_testchannel->GetY()[i-1]));
    }
    g_offset_correction->SetMarkerStyle(20);
    g_offset_correction->SetMarkerSize(0.5);
    g_offset_correction->SetMarkerColor(1);
    g_offset_correction->SetLineColor(1);
    g_offset_correction->SetTitle("Offset correction;iteration;Offset correction [ns]");
    g_offset_correction->Draw("APL");
    TF1* series = new TF1("series", "1/(x^2)", 1, maxIteration);
    series->SetLineColor(kRed);
    series->Draw("same");
    TGraph* g_offset_correction_ratio = new TGraph();
    for(int i = 1; i < g_offset_correction->GetN(); i++) {
        g_offset_correction_ratio->SetPoint(g_offset_correction_ratio->GetN(), i, g_offset_correction->GetY()[i]/g_offset_correction->GetY()[i-1]);
    }
    g_offset_correction_ratio->SetMarkerStyle(20);
    g_offset_correction_ratio->SetMarkerSize(0.5);
    g_offset_correction_ratio->SetMarkerColor(1);
    g_offset_correction_ratio->SetLineColor(1);
    g_offset_correction_ratio->SetTitle("Offset correction ratio;iteration;Offset correction ratio");
    g_offset_correction_ratio->Draw("APL");

    if (saveImage)
        c_offset->SaveAs(Form("TimeOffsetCalibEvolution_%s.png", reducedFilename.c_str()));

    TCanvas* c_offsetHist = new TCanvas("c_offsetHist", "c_offsetHist", 1800, 500);
    c_offsetHist->Divide(2,1);
    c_offsetHist->cd(1);
    h_inputOffsets->SetLineColor(kRed);
    h_inputOffsets->SetLineStyle(2);
    h_inputOffsets->SetLineWidth(4);
    h_outputOffsets->SetLineColor(kBlue);
    h_inputOffsets->Draw();
    h_outputOffsets->Draw("same");
    TLegend *inputputputleg = new TLegend(0.1,0.7,0.48,0.9);
    inputputputleg->AddEntry(h_inputOffsets, "Input offsets", "l");
    inputputputleg->AddEntry(h_outputOffsets, "Output offsets", "l");
    inputputputleg->Draw();
    c_offsetHist->cd(2);
    g_inputOffsets->SetMarkerStyle(20);
    g_inputOffsets->SetMarkerSize(0.5);
    g_inputOffsets->SetMarkerColor(kRed);
    g_inputOffsets->SetLineColor(kRed);
    g_outputOffsets->SetMarkerStyle(20);
    g_outputOffsets->SetMarkerSize(0.5);
    g_outputOffsets->SetMarkerColor(kBlue);
    g_outputOffsets->SetLineColor(kBlue);
    g_inputOffsets->SetTitle("Input vs. output offsets;Hits in channel;Offset [ns]");
//  g_inputOffsets->Draw("APl");
//  g_outputOffsets->Draw("Pl");
    g_residuals->SetMarkerStyle(20);
    g_residuals->SetMarkerSize(0.5);
    g_residuals->SetMarkerColor(kBlack);
    g_residuals->SetTitle("Residuals;Channel;Offset_{output} - Offset_{input} [ps]");
    g_residuals->Draw("APL");

    if(saveImage)
        c_offsetHist->SaveAs(Form("TimeOffsetCalibInputOutput_%s.png", reducedFilename.c_str()));

    TCanvas* c_offsetCorrection = new TCanvas("c_offsetCorrection", "c_offsetCorrection", 1800, 600);
    h_offset_correction_ratio_last->Draw();

    if(saveImage)
        c_offsetCorrection->SaveAs(Form("TimeOffsetCalibOffsetCorrection_%s.png", reducedFilename.c_str()));

    TCanvas* c_delta12 = new TCanvas("c_delta12", "c_delta12", 1800, 600);
    c_delta12->Divide(2,1);
    c_delta12->cd(1);
    g_magnitudeDeltaTk->SetMarkerStyle(20);
    g_magnitudeDeltaTk->SetMarkerSize(0.5);
    g_magnitudeDeltaTk->SetMarkerColor(1);
    g_magnitudeDeltaTk->SetLineColor(1);
    g_magnitudeDeltaTk->SetTitle("|#Delta T^{k}|; k; |#Delta T^{k}|");
    g_magnitudeDeltaTk->Draw("APL");
//  c_delta12->cd(2);
//  g_delta_t2->SetMarkerStyle(20);
//  g_delta_t2->SetMarkerSize(0.5);
//  g_delta_t2->SetMarkerColor(2);
//  g_delta_t2->SetLineColor(2);
//  g_delta_t2->SetTitle("#sum_n (t_2-s_2);Iteration;Delta t2 [ns]");
//  g_delta_t2->Draw("APL");
//  c_delta12->cd(1);
//  g_delta_t2->Draw("samePL");
//  TGraph* g_deltat1_sum = new TGraph();
//  double sum = 0;
//  for(int i = 0; i < g_delta_t1->GetN(); i++) {
//    sum += g_delta_t1->GetY()[i] - g_delta_t2->GetY()[i];
//    g_deltat1_sum->SetPoint(g_deltat1_sum->GetN(), i, sum);
//  }
//  g_deltat1_sum->SetMarkerStyle(20);
//  g_deltat1_sum->SetMarkerSize(0.5);
//  g_deltat1_sum->SetMarkerColor(3);
//  g_deltat1_sum->SetLineColor(3);
//  g_deltat1_sum->Draw("samePL");
//  TGraph* g_computed_offset_testchannel_new = new TGraph(*g_computed_offset_testchannel);
//  g_computed_offset_testchannel_new->Draw("samePL");
//  g_computed_offset_testchannel_new->SetMarkerColor(4);
//  g_computed_offset_testchannel_new->SetLineColor(4);
//  line->Draw("same");
//  g_delta_t1t2ratio->SetMarkerStyle(20);
//  g_delta_t1t2ratio->SetMarkerSize(0.5);
//  g_delta_t1t2ratio->SetMarkerColor(kMagenta);
//  g_delta_t1t2ratio->SetLineColor(kMagenta);
//  g_delta_t1t2ratio->Draw("APL");
//  g_delta_t1t2ratio->SetTitle("1/N_{A}#sum_{n} #frac{t_{2}-s_{2}}{t_{1}-s_{1}};Iteration;Delta t2/Delta t1");
    c_delta12->cd(2);
    g_magnitudeDeltaTk_ratio->SetMarkerStyle(20);
    g_magnitudeDeltaTk_ratio->SetMarkerSize(0.5);
    g_magnitudeDeltaTk_ratio->SetMarkerColor(kMagenta);
    g_magnitudeDeltaTk_ratio->SetLineColor(kMagenta);
    g_magnitudeDeltaTk_ratio->Draw("APL");
    g_magnitudeDeltaTk_ratio->SetTitle("|#Delta T^{k}|/|#Delta T^{k-1}|; k; |#Delta T^{k}|/|#Delta T^{k-1}|");

    if(saveImage)
        c_delta12->SaveAs(Form("TimeOffsetCalibDelta12_%s.png", reducedFilename.c_str()));


    if (batchMode) {
        gROOT->SetBatch(kFALSE);
    }
    return sigma_offset;

}

