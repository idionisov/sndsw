#pragma once

#include "sndSciFiBaseCut.h"
#include "TChain.h"

namespace snd {
  namespace analysis_cuts {
    class maxSciFiHitsCut : public snd::analysis_cuts::sciFiBaseCut {
    private :
      int hitThreshold;
    public :
      maxSciFiHitsCut(int threshold, TChain * ch);
      ~maxSciFiHitsCut(){;}

      bool passCut();
    };
  }
}
