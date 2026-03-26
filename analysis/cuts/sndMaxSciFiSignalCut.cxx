#include "sndMaxSciFiSignalCut.h"
#include <numeric>

namespace snd::analysis_cuts {

  maxSciFiSignalCut::maxSciFiSignalCut(double threshold, TChain * ch) : sciFiBaseCut(ch), signalThreshold(threshold) {
    cutName = "Maximum Total SciFi Signal Cut";
    shortName = "maxSciFiSignal";
    nbins.push_back(1);
    range_start.push_back(0);
    range_end.push_back(1);
  }

  bool maxSciFiSignalCut::passCut() {
    initializeEvent();

    double totalSignalV = std::accumulate(signal_per_plane_vertical.begin(), signal_per_plane_vertical.end(), 0.);
    double totalSignalH = std::accumulate(signal_per_plane_horizontal.begin(), signal_per_plane_horizontal.end(), 0.);

    return (totalSignalV + totalSignalH < signalThreshold);
  }
}
