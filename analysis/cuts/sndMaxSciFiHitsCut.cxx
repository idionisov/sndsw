#include "sndMaxSciFiHitsCut.h"
#include <numeric>

namespace snd::analysis_cuts {

  maxSciFiHitsCut::maxSciFiHitsCut(int threshold, TChain * ch) : sciFiBaseCut(ch), hitThreshold(threshold) {
    cutName = "Maximum Total SciFi Hits Cut";
    shortName = "maxSciFiHits";
    nbins.push_back(1);
    range_start.push_back(0);
    range_end.push_back(1);
  }

  bool maxSciFiHitsCut::passCut() {
    initializeEvent();

    int nHitsV = std::accumulate(hits_per_plane_vertical.begin(), hits_per_plane_vertical.end(), 0);
    int nHitsH = std::accumulate(hits_per_plane_horizontal.begin(), hits_per_plane_horizontal.end(), 0);

    return (nHitsV + nHitsH < hitThreshold);
  }
}
