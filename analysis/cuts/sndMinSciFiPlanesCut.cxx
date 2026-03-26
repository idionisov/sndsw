#include "sndMinSciFiPlanesCut.h"
#include <algorithm>

namespace snd::analysis_cuts {

  minSciFiPlanesCut::minSciFiPlanesCut(int minH, int minV, TChain * ch) : sciFiBaseCut(ch), minPlanesH(minH), minPlanesV(minV) {
    cutName = "Minimum SciFi Planes Cut";
    shortName = "minSciFiPlanes";
    nbins.push_back(1);
    range_start.push_back(0);
    range_end.push_back(1);
  }

  bool minSciFiPlanesCut::passCut() {
    initializeEvent();

    int nPlanesV = std::count_if(hits_per_plane_vertical.begin(), hits_per_plane_vertical.end(), [](int n){ return n > 0; });
    int nPlanesH = std::count_if(hits_per_plane_horizontal.begin(), hits_per_plane_horizontal.end(), [](int n){ return n > 0; });

    return (nPlanesV >= minPlanesV && nPlanesH >= minPlanesH);
  }
}
