// +build ignore

// Std libs
#include <stdio.h>

// ROOT libs
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

// topreco lib
#include "topreco.h"


void test_root() {
	//gSystem->Load("topreco");

	auto f = TFile::Open("../testdata/data.root");
	auto t = (TTree*)f->Get("nominal");

	// Init the TopReconstruction tool
	auto smearFile = "../testdata/smearingHistos.root";
	auto nSmearing = 100;
	auto debugMode = false;
	initTopReconstruction(smearFile, nSmearing, debugMode);
	
	const int NLEPS = 10;
	const int NJETS = 10;

	Long64_t evtNum;
	Int_t nlep;
	float lepPt[NLEPS];
	float lepEta[NLEPS];
	float lepPhi[NLEPS];
	Int_t lepPID[NLEPS];

	Int_t njet;
	Int_t nbjet;
	float jetE[NJETS];
	float jetPt[NJETS];
	float jetEta[NJETS];
	float jetPhi[NJETS];
	float jetMV2[NJETS];
	float metMet, metPhi;

	t->SetBranchAddress("eventNumber", &evtNum);
	t->SetBranchAddress("d_nlep", &nlep);
	t->SetBranchAddress("d_lep_pt", &lepPt);
	t->SetBranchAddress("d_lep_eta", &lepEta);
	t->SetBranchAddress("d_lep_phi", &lepPhi);
	t->SetBranchAddress("d_lep_pid", &lepPID);

	t->SetBranchAddress("d_njet", &njet);
	t->SetBranchAddress("d_nbjet", &nbjet);
	t->SetBranchAddress("d_jet_e", &jetE);
	t->SetBranchAddress("d_jet_pt", &jetPt);
	t->SetBranchAddress("d_jet_eta", &jetEta);
	t->SetBranchAddress("d_jet_phi", &jetPhi);
	t->SetBranchAddress("d_jet_mv2c10", &jetMV2);

	t->SetBranchAddress("d_met_met", &metMet);
	t->SetBranchAddress("d_met_phi", &metPhi);

	Int_t N = t->GetEntries();
	printf("entries: %d\n", N);

	const bool debug = false;

	Int_t nok = 0;

	for (Int_t n = 0; n < N; n++) {
		t->GetEntry(n);
		if (n%10 == 0) {
			printf("entry: %d\n", n);
		}
		if (debug) {
			printf("nlep: %d\n", nlep);
			for (int i = 0; i < nlep; i++) {
				printf("lep[%d]: (%g, %g, %g, %d)\n",
						i, lepPt[i], lepEta[i], lepPhi[i], lepPID[i]
				);
			}
			printf("njet: %d\n", njet);
			for (int i = 0; i < njet; i++) {
				printf("jet[%d]: (%g, %g, %g, %g)\n",
						i, jetPt[i], jetEta[i], jetPhi[i], jetMV2[i]
				);
			}
		}

		double etx = double(metMet)*TMath::Cos(metPhi);
		double ety = double(metMet)*TMath::Sin(metPhi);

		Int_t j0b = jetMV2[0] > 0.691 ? 1 : 0;
		Int_t j1b = jetMV2[1] > 0.691 ? 1 : 0;

		TLorentzVector lep0, lep1, jet0, jet1, t, tbar;
		lep0.SetPtEtaPhiM(lepPt[0], lepEta[0], lepPhi[0], 0);
		lep1.SetPtEtaPhiM(lepPt[1], lepEta[1], lepPhi[1], 0);
		jet0.SetPtEtaPhiM(jetPt[0], jetEta[0], jetPhi[0], jetE[0]);
		jet1.SetPtEtaPhiE(jetPt[1], jetEta[1], jetPhi[1], jetE[1]);
		Int_t pdg0 = lepPID[0];
		Int_t pdg1 = lepPID[1];
		if (lepPID[0] > 0) {
			std::swap(lep0, lep1);
			std::swap(pdg0, pdg1);
		}

		// Call the reconstruction
		Int_t nIter = runTopReconstruction(lep0, lep1, lepPID[0], lepPID[1],
						   jet0, jet1, j0b, j1b,
						   double(etx), double(ety),
						   &t, &tbar);
		
		if (nIter > 0) {
		  printf("evt: %lld (entry=%d)\n"  , evtNum, n);
		  printf("Iterations: %d\n", nIter);
		  printf("t         : (%g, %g, %g, %g)\n", t.Px(), t.Py(), t.Pz(), t.E());
		  printf("tbar      : (%g, %g, %g, %g)\n", tbar.Px(), tbar.Py(), tbar.Pz(), tbar.E());
		  nok++;
		}
	}

	printf("n-reco: %d\n", nok);
}
