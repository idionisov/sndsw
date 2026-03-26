#include "sndHasVetoHitsCut.h"

#include "TClonesArray.h"
#include "TChain.h"
#include "MuFilterHit.h"

namespace snd::analysis_cuts {

  hasVetoHitsCut::hasVetoHitsCut(TChain * ch) : MuFilterBaseCut(ch) {
    cutName = "At least one hit in veto";

    shortName = "HasVetoHits";
    nbins = std::vector<int>{16};
    range_start = std::vector<double>{0};
    range_end = std::vector<double>{16};
    plot_var = std::vector<double>{-1};

  }

  bool hasVetoHitsCut::passCut(){
    MuFilterHit * hit;
    TIter hitIterator(muFilterDigiHitCollection);
    
    plot_var[0] = 0;


    while ( (hit = (MuFilterHit*) hitIterator.Next()) ){
      if (hit->GetSystem() == 1) plot_var[0] += 1;
    }

    if (plot_var[0] > 0) return true;
    return false;
  }
}
