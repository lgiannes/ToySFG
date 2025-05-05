

namespace common_utils{
    typedef vector <vector<double>> Matrix;

    inline void draw_matrix(Matrix C,int nChannels, string title) {

      TH2D *hC = new TH2D("hC", ";Channel i;Channel j",
                          nChannels, 0, nChannels,
                          nChannels, 0, nChannels);

      for (int i = 0; i < nChannels; ++i) {
        for (int j = 0; j < nChannels; ++j) {
          hC->SetBinContent(i + 1, j + 1, C[i][j]);  // Bin numbers start at 1
        }
      }

      hC->SetStats(0);
      hC->Draw("COLZ");
      hC->SetTitle(title.c_str());

      return;
    }
};