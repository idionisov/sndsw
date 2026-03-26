#include "sndTridentDensityCut.h"
#include "sndSciFiTools.h"
#include "TClonesArray.h"
#include "sndScifiHit.h"

namespace snd::analysis_cuts {

  tridentDensityCut::tridentDensityCut(int r, int t, TChain * ch) : sciFiBaseCut(ch), radius(r), threshold(t) {
    cutName = "Trident SciFi Density Cut";
    shortName = "tridentDensity";
    nbins.push_back(1);
    range_start.push_back(0);
    range_end.push_back(1);
  }

  bool tridentDensityCut::passCut() {
    initializeEvent();

    double total_density = 0;
    sndScifiHit * hit;
    TIter hitIterator(scifiDigiHitCollection);

    while ( (hit = (sndScifiHit*) hitIterator.Next()) ){
      if (hit->isValid()){
        total_density += snd::analysis_tools::densityScifi(hit->GetChannelID(), *scifiDigiHitCollection, radius, 1000000, false);
      }
    }

    return (total_density < threshold);
  }
}
