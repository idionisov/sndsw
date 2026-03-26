#pragma once

#include "sndMuFilterBaseCut.h"

#include "TChain.h"

namespace snd {
  namespace analysis_cuts {
  
    class hasVetoHitsCut : public snd::analysis_cuts::MuFilterBaseCut {
    public :
      hasVetoHitsCut(TChain * ch);
      ~hasVetoHitsCut(){;}
      bool passCut();
    };

  }
}
