#pragma once

#include "sndSciFiBaseCut.h"
#include "TChain.h"

namespace snd {
  namespace analysis_cuts {
    class maxPlaneSciFiHitsCut : public snd::analysis_cuts::sciFiBaseCut {
    private :
      int planeThreshold;
    public :
      maxPlaneSciFiHitsCut(int threshold, TChain * ch);
      ~maxPlaneSciFiHitsCut(){;}

      bool passCut();
    };
  }
}
