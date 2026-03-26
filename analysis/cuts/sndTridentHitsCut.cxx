#include "sndTridentHitsCut.h"
#include <numeric>

namespace snd::analysis_cuts {

  tridentHitsCut::tridentHitsCut(int threshold1, int threshold2, TChain * ch) : sciFiBaseCut(ch), hitThreshold1(threshold1), hitThreshold2(threshold2) {
    cutName = "Trident SciFi Hits Cut";
    shortName = "tridentSciFiHits";
    nbins.push_back(1);
    range_start.push_back(0);
    range_end.push_back(1);
  }

  bool tridentHitsCut::passCut() {
    initializeEvent();

    int nHitsV = std::accumulate(hits_per_plane_vertical.begin(), hits_per_plane_vertical.end(), 0);
    int nHitsH = std::accumulate(hits_per_plane_horizontal.begin(), hits_per_plane_horizontal.end(), 0);

    return (nHitsV >= hitThreshold1 && nHitsH >= hitThreshold2) || (nHitsH >= hitThreshold1 && nHitsV >= hitThreshold2);
  }
}
