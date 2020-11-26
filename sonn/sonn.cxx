#include "sonn.h"
#include <vector>

#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TMath.h"
#include "Math/Polynomial.h"
#include "TH1F.h"

namespace AnalyticalTopReconstruction{

  std::vector<TLorentzVector> Sonnenschein(TLorentzVector tlv_lep, TLorentzVector tlv_lepbar, int pdgId_lep, int pdgId_lepbar, TLorentzVector tlv_jet, TLorentzVector tlv_jetbar, Double_t Emiss_x, Double_t Emiss_y, int nBJets);

}

tlv_t get_tlv(tlv_t *tlvs, int i) {
	return tlvs[i];
}

tlvs_t sonn(
		tlv_t lep, tlv_t lepbar, int pdgID_lep, int pdgID_lepbar,
		tlv_t jet, tlv_t jetbar, 
		double emissx, double emissy,
		int nbjets) {

	TLorentzVector tlvLep;
	tlvLep.SetPxPyPzE(lep.Px, lep.Py, lep.Pz, lep.E);
	TLorentzVector tlvLepBar;
	tlvLepBar.SetPxPyPzE(lepbar.Px, lepbar.Py, lepbar.Pz, lepbar.E);

	TLorentzVector tlvJet;
	tlvJet.SetPxPyPzE(jet.Px, jet.Py, jet.Pz, jet.E);
	TLorentzVector tlvJetBar;
	tlvJetBar.SetPxPyPzE(jetbar.Px, jetbar.Py, jetbar.Pz, jetbar.E);

	auto out = AnalyticalTopReconstruction::Sonnenschein(
			tlvLep, tlvLepBar, pdgID_lep, pdgID_lepbar,
			tlvJet, tlvJetBar, emissx, emissy, nbjets);

	tlvs_t tlvs;
	tlvs.n = out.size();
	tlvs.tlvs = (tlv_t*)malloc(sizeof(tlv_t)*tlvs.n);
	for (int i = 0; i < tlvs.n; i++) {
		tlvs.tlvs[i].Px = out[i].Px();
		tlvs.tlvs[i].Py = out[i].Py();
		tlvs.tlvs[i].Pz = out[i].Pz();
		tlvs.tlvs[i].E = out[i].E();
	}

	return tlvs;
}

