extern "C" {
#include "./libtopreco.h"
}
R__LOAD_LIBRARY(./libtopreco.so)

#include "TLorentzVector.h"

P4_t newP4(TLorentzVector tlv){
	return P4_t{tlv.Px(), tlv.Py(), tlv.Pz(), tlv.E()};
}

Int_t runTopReconstruction(TLorentzVector l, TLorentzVector lbar, Int_t pid, Int_t pidbar,
			   TLorentzVector j, TLorentzVector jbar, Int_t jbtag, Int_t jbarbtag,
			   double etx, double ety,
			   TLorentzVector *t, TLorentzVector *tbar) {  
  P4_t top0, top1;
  Int_t ok = ReconstructTops(
			     newP4(l), newP4(lbar), pid, pidbar,
			     newP4(j), newP4(jbar), jbtag, jbarbtag,
			     etx, ety, &top0, &top1);
  
  t->SetPxPyPzE(top0.Px, top0.Py, top0.Pz, top0.E);
  tbar->SetPxPyPzE(top1.Px, top1.Py, top1.Pz, top1.E);
    
  return ok;
}
