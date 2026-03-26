#pragma once

#include "sndSciFiBaseCut.h"
#include "TChain.h"

namespace snd {
  namespace analysis_cuts {
    class tridentDensityCut : public snd::analysis_cuts::sciFiBaseCut {
    private :
      int radius;
      int threshold;
    public :
      tridentDensityCut(int r, int t, TChain * ch);
      ~tridentDensityCut(){;}

      bool passCut();
    };
  }
}
