extern "C" {
#include "./libtopreco.h"
}
R__LOAD_LIBRARY(./libtopreco.so)

#include "TLorentzVector.h"

// newP4 converts TLorentzVector into a type understood by the tool.
P4_t newP4(TLorentzVector tlv){
	return P4_t{tlv.Px(), tlv.Py(), tlv.Pz(), tlv.E()};
}

// initTopReconstruction initializes the top reconstruction tool
void initTopReconstruction(string smearFile, int nSmearing, bool debug) {
  InitTopBuilder((char*)smearFile.c_str(), nSmearing, int(debug));
}

// runTopReconstruction performs the reconstruction and fill (t, tbar) 4-momentum
// This function returns the number of successful smearing iteration.
Int_t runTopReconstruction(TLorentzVector l, TLorentzVector lbar, Int_t pid, Int_t pidbar,
			   TLorentzVector j, TLorentzVector jbar, Int_t jbtag, Int_t jbarbtag,
			   double etx, double ety, TLorentzVector *t, TLorentzVector *tbar) {  
  P4_t top0, top1;
  Int_t nIter = ReconstructTops(
			     newP4(l), newP4(lbar), pid, pidbar,
			     newP4(j), newP4(jbar), jbtag, jbarbtag,
			     etx, ety, &top0, &top1);
  
  t->SetPxPyPzE(top0.Px, top0.Py, top0.Pz, top0.E);
  tbar->SetPxPyPzE(top1.Px, top1.Py, top1.Pz, top1.E);
    
  return nIter;
}
