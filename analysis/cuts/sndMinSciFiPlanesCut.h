#pragma once

#include "sndSciFiBaseCut.h"
#include "TChain.h"

namespace snd {
  namespace analysis_cuts {
    class minSciFiPlanesCut : public snd::analysis_cuts::sciFiBaseCut {
    private :
      int minPlanesH;
      int minPlanesV;
    public :
      minSciFiPlanesCut(int minH, int minV, TChain * ch);
      ~minSciFiPlanesCut(){;}

      bool passCut();
    };
  }
}
