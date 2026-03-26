#pragma once

#include "sndSciFiBaseCut.h"
#include "TChain.h"

namespace snd {
  namespace analysis_cuts {
    class maxPlaneSciFiSignalCut : public snd::analysis_cuts::sciFiBaseCut {
    private :
      double planeThreshold;
    public :
      maxPlaneSciFiSignalCut(double threshold, TChain * ch);
      ~maxPlaneSciFiSignalCut(){;}

      bool passCut();
    };
  }
}
