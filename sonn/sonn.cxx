#include "sonn.h"
#include <vector>

#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TMath.h"
#include "Math/Polynomial.h"
#include "TH1F.h"
#include "TFile.h"

// Smearing histograms
extern TH1F *h_pT_lep_ee_smear;
extern TH1F *h_pT_lep_mm_smear;
extern TH1F *h_mlblb_smear;
extern TH1F *h_theta_jet_smear;
extern TH1F *h_theta_jetbar_smear;
extern TH1F *h_pT_jet_smear;
extern TH1F *h_theta_lep_ee_smear;
extern TH1F *h_theta_lep_mm_smear;

void load_smearing_histos(const char *fname) {
	auto f = TFile::Open(fname);
	h_pT_lep_ee_smear =    (TH1F*)f->Get("h_pT_lep_ee_smear");
	h_pT_lep_mm_smear =    (TH1F*)f->Get("h_pT_lep_mm_smear");
	h_mlblb_smear =        (TH1F*)f->Get("mlblb_smear");
	h_theta_jet_smear =    (TH1F*)f->Get("h_theta_jet_smear");
	h_theta_jetbar_smear = (TH1F*)f->Get("h_theta_jetbar_smear");
	h_pT_jet_smear =       (TH1F*)f->Get("h_pT_jet_smear");
	h_theta_lep_ee_smear = (TH1F*)f->Get("h_theta_lep_ee_smear");
	h_theta_lep_mm_smear = (TH1F*)f->Get("h_theta_lep_mm_smear");
}

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

