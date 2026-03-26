#include "sndMaxPlaneSciFiSignalCut.h"
#include <algorithm>

namespace snd::analysis_cuts {

  maxPlaneSciFiSignalCut::maxPlaneSciFiSignalCut(double threshold, TChain * ch) : sciFiBaseCut(ch), planeThreshold(threshold) {
    cutName = "Maximum Plane SciFi Signal Cut";
    shortName = "maxPlaneSciFiSignal";
    nbins.push_back(1);
    range_start.push_back(0);
    range_end.push_back(1);
  }

  bool maxPlaneSciFiSignalCut::passCut() {
    initializeEvent();

    for (double s : signal_per_plane_vertical) {
      if (s >= planeThreshold) return false;
    }
    for (double s : signal_per_plane_horizontal) {
      if (s >= planeThreshold) return false;
    }

    return true;
  }
}
