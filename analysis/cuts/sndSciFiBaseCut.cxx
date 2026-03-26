#include "sndSciFiBaseCut.h"

#include <vector>

#include "TClonesArray.h"
#include "TChain.h"
#include "sndScifiHit.h"

namespace snd::analysis_cuts {

  TClonesArray * sciFiBaseCut::scifiDigiHitCollection = 0;
  TChain * sciFiBaseCut::tree = 0;
  unsigned long int sciFiBaseCut::read_entry = -1;

  std::vector<int> sciFiBaseCut::hits_per_plane_vertical = std::vector<int>(5, 0);
  std::vector<int> sciFiBaseCut::hits_per_plane_horizontal = std::vector<int>(5, 0);
  std::vector<double> sciFiBaseCut::signal_per_plane_vertical = std::vector<double>(5, 0.);
  std::vector<double> sciFiBaseCut::signal_per_plane_horizontal = std::vector<double>(5, 0.);

  sciFiBaseCut::sciFiBaseCut(TChain * ch){
    if (tree == 0){
      tree = ch;
      scifiDigiHitCollection = new TClonesArray("sndScifiHit", 3000);
      tree->SetBranchAddress("Digi_ScifiHits", &scifiDigiHitCollection);
    }
  }

  void sciFiBaseCut::initializeEvent(){
    if (read_entry != tree->GetReadEntry()){
      read_entry = tree->GetReadEntry();

      // Clear hits and signal vectors
      std::fill(hits_per_plane_vertical.begin(), hits_per_plane_vertical.end(), 0);
      std::fill(hits_per_plane_horizontal.begin(), hits_per_plane_horizontal.end(), 0);
      std::fill(signal_per_plane_vertical.begin(), signal_per_plane_vertical.end(), 0.);
      std::fill(signal_per_plane_horizontal.begin(), signal_per_plane_horizontal.end(), 0.);

      // Add valid hits to vectors
      sndScifiHit * hit;
      TIter hitIterator(scifiDigiHitCollection);

      while ( (hit = (sndScifiHit*) hitIterator.Next()) ){
	if (hit->isValid()){
	  int sta = hit->GetStation();
	  if (hit->isVertical()){
	    hits_per_plane_vertical[sta-1]++;
	    signal_per_plane_vertical[sta-1] += hit->GetSignal();
	  } else {
	    hits_per_plane_horizontal[sta-1]++;
	    signal_per_plane_horizontal[sta-1] += hit->GetSignal();
	  }
	}
      }
    }
  }
}
