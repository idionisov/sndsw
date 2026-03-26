#pragma once

#include "sndSciFiBaseCut.h"
#include "TChain.h"

namespace snd {
  namespace analysis_cuts {
    class maxSciFiSignalCut : public snd::analysis_cuts::sciFiBaseCut {
    private :
      double signalThreshold;
    public :
      maxSciFiSignalCut(double threshold, TChain * ch);
      ~maxSciFiSignalCut(){;}

      bool passCut();
    };
  }
}
