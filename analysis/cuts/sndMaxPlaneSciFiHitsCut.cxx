#include "sndMaxPlaneSciFiHitsCut.h"
#include <algorithm>

namespace snd::analysis_cuts {

  maxPlaneSciFiHitsCut::maxPlaneSciFiHitsCut(int threshold, TChain * ch) : sciFiBaseCut(ch), planeThreshold(threshold) {
    cutName = "Maximum Plane SciFi Hits Cut";
    shortName = "maxPlaneSciFiHits";
    nbins.push_back(1);
    range_start.push_back(0);
    range_end.push_back(1);
  }

  bool maxPlaneSciFiHitsCut::passCut() {
    initializeEvent();

    for (int n : hits_per_plane_vertical) {
      if (n >= planeThreshold) return false;
    }
    for (int n : hits_per_plane_horizontal) {
      if (n >= planeThreshold) return false;
    }

    return true;
  }
}
