#include <stdio.h> 
#include <iostream>

#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "TMath.h"
#include "Math/Polynomial.h"
#include "TH1F.h"
#include "TFile.h"

// Smearing histograms
TH1F *h_pT_lep_ee_smear;
TH1F *h_pT_lep_mm_smear;
TH1F *h_mlblb_smear;
TH1F *h_theta_jet_smear;
TH1F *h_theta_jetbar_smear;
TH1F *h_pT_jet_smear;
TH1F *h_theta_lep_ee_smear;
TH1F *h_theta_lep_mm_smear;

namespace AnalyticalTopReconstruction{
  
  // Compute the complete four vectors of top quarks
  std::vector<TLorentzVector> Sonnenschein(TLorentzVector tlv_lep, TLorentzVector tlv_lepbar, int pdgId_lep, int pdgId_lepbar,
					   TLorentzVector tlv_jet, TLorentzVector tlv_jetbar, bool isb_jet, bool isb_jetbar, 
					   Double_t Emiss_x, Double_t Emiss_y) {
    
    // Debug
    bool print = true;
    
    // Smearing parameter
    int Nsmear = 10;
    bool applySmearing = true;
    if (!applySmearing) Nsmear=1; 
    
    // Container for output 4-momentum
    std::vector<TLorentzVector> return_top_tbar;
    return_top_tbar.push_back(TLorentzVector());
    return_top_tbar.push_back(TLorentzVector());

    // Run the reconstruction only if the two jets are btags
    if (!isb_jet || !isb_jetbar)
      return return_top_tbar;
       
    bool lep_is_e     = ( abs(pdgId_lep)    == 11 );
    bool lep_is_mu    = ( abs(pdgId_lep)    == 13 );
    bool lepbar_is_e  = ( abs(pdgId_lepbar) == 11 );
    bool lepbar_is_mu = ( abs(pdgId_lepbar) == 13 );
    
    // the lepton TLV is still in MeV, rescale!
    // tlv_lep.SetPtEtaPhiE(   tlv_lep.Pt()/1000.,    tlv_lep.Eta(),    tlv_lep.Phi(),    tlv_lep.E()/1000.);
    // tlv_lepbar.SetPtEtaPhiE(tlv_lepbar.Pt()/1000., tlv_lepbar.Eta(), tlv_lepbar.Phi(), tlv_lepbar.E()/1000.);
    
    // particle masses (in GeV)
    Double_t m_W      = 80.379;
    Double_t m_Wbar   = 80.379;
    Double_t m_b      = 4.18;
    Double_t m_bbar   = 4.18;
    Double_t m_nu     = 0.;
    Double_t m_nubar  = 0.;
    Double_t m_top    = 172.5;
    Double_t m_topbar = 172.5;
    Double_t m_lep    = 0;
    Double_t m_lepbar = 0;
    
    TLorentzVector nu, nubar;
    TLorentzVector nu_NW, nubar_NW;
    TLorentzVector jet, jetbar;
    TLorentzVector lep,lepbar;
    TLorentzVector lep_pt_smear,lepbar_pt_smear;
    TLorentzVector lep_nosmear,lepbar_nosmear;
    TLorentzVector jet_pt_smear,jetbar_pt_smear;
    TLorentzVector jet_nosmear,jetbar_nosmear;
    TLorentzVector top,topbar;
    TLorentzVector Top_calc, Topbar_calc;
    TLorentzVector nu_calc, nubar_calc;
    TLorentzVector lplus;
    TLorentzVector lminus;
    TLorentzVector lep0, lepbar0;
    TLorentzVector lep1, lepbar1;
    TLorentzVector Top1,Topbar1,Top2;

    
    Int_t i_analyzed_event; // Not really used
    Int_t n_solutions;
    
    Double_t p_nu_y = -999.;
    Double_t p_nu_z = -999.;
    Double_t p_nubar_x = -999.;
    Double_t p_nubar_y = -999.;
    Double_t p_nubar_z = -999.;
    Double_t delta_z = 100000.;
    Double_t p_nu_x_close = 100000.;
    Double_t p_nu_x_close_f = -999.;
    Double_t p_nu_y_f = -999.;
    Double_t p_nu_z_f = -999.;
    Double_t p_nubar_x_f = -999.;
    Double_t p_nubar_y_f = -999.;
    Double_t p_nubar_z_f = -999.;
    Double_t delta_TOP_mass = 100000.;
    Double_t Top_mass_f = 100000.;
    Double_t Topbar_mass_f = 100000.;
    Double_t Top_mass_f_1 = 100000.;
    std::vector<Double_t> mtt_val;
	
    Double_t Vec_Top[2][3];
    Double_t Vec_Topbar[2][3];
    std::vector<Double_t> weights_com;
    for (int i_n1=0; i_n1<2; i_n1++ ) {
      for (int i_n2=0; i_n2<3; i_n2++ ) {
	Vec_Top[i_n1][i_n2]    = 0.;
	Vec_Topbar[i_n1][i_n2] = 0.;
      }
    }
    
    // loop over jets
    for (int i_jets=0; i_jets<2; i_jets++) {
      
      Double_t weight_s_sum = 0.;
      TRandom3 r;

      if (print) {
	std::cout << " Jets combination " << i_jets << ": " << std::endl;
      }
      
      int nIterations = 0;
      for (int i_smear=0; i_smear<Nsmear; i_smear++) {
	
	Double_t smear_scale_0     = 1.;
	Double_t smear_scale_1     = 1.;
	Double_t smear_scale_jet_0 = 1.;
	Double_t smear_scale_jet_1 = 1.;

	if (applySmearing) {
	  if (lep_is_e)          smear_scale_0 = h_pT_lep_ee_smear->GetRandom();
	  else if (lep_is_mu)    smear_scale_0 = h_pT_lep_mm_smear->GetRandom();
	  else                   smear_scale_0 = 1.;
	  if (lepbar_is_e)       smear_scale_1 = h_pT_lep_ee_smear->GetRandom();
	  else if (lepbar_is_mu) smear_scale_1 = h_pT_lep_mm_smear->GetRandom();
	  else                   smear_scale_1 = 1.;
	}
	
	lepbar_pt_smear = tlv_lepbar;
	lepbar_pt_smear.SetPtEtaPhiM(lepbar_pt_smear.Pt()*smear_scale_0, lepbar_pt_smear.Eta(), lepbar_pt_smear.Phi(), m_lepbar);
	lep_pt_smear    = tlv_lep;
	lep_pt_smear.SetPtEtaPhiM(lep_pt_smear.Pt()*smear_scale_1,lep_pt_smear.Eta(),lep_pt_smear.Phi(),m_lep);
	lepbar_nosmear  = tlv_lepbar;
	lep_nosmear     = tlv_lep;
	
	// let's define axis transverse to lepton momentun (for theta smearing accounting)
	TVector3 lep_v_pt_1    = lep_pt_smear.Vect();
	TVector3 lepbar_v_pt_1 = lepbar_pt_smear.Vect();
	
	TVector3 beam_axis1;
	beam_axis1.SetXYZ(0,0,1);
	
	TVector3 lep_axis = lep_v_pt_1;
	TVector3 transe_lep_axis = beam_axis1.Cross(lep_axis);
	
	TVector3 lepbar_axis = lepbar_v_pt_1;
	TVector3 transe_lepbar_axis = beam_axis1.Cross(lepbar_axis);
	
	TVector3 lep_v_pt = lep_v_pt_1;
	TVector3 lepbar_v_pt = lepbar_v_pt_1;
	
	Double_t smear_angle_ee_lep    = 0.;
	Double_t smear_angle_ee_lepbar = 0.;
	Double_t smear_angle_mm_lep    = 0.;
	Double_t smear_angle_mm_lepbar = 0.;

	if (applySmearing) {
	  smear_angle_ee_lep    = h_theta_lep_ee_smear->GetRandom();
	  smear_angle_ee_lepbar = h_theta_lep_ee_smear->GetRandom();
	  smear_angle_mm_lep    = h_theta_lep_mm_smear->GetRandom();
	  smear_angle_mm_lepbar = h_theta_lep_mm_smear->GetRandom();
	}
	
	if (tlv_lep.Pt() > tlv_lepbar.Pt()) {
	  if (lep_is_e)          lep_v_pt.Rotate(smear_angle_ee_lep, transe_lep_axis);
	  else if (lep_is_mu)    lep_v_pt.Rotate(smear_angle_mm_lep, transe_lep_axis);
	}
	else {
	  if (lepbar_is_e)       lepbar_v_pt.Rotate(smear_angle_ee_lepbar, transe_lepbar_axis);
	  else if (lepbar_is_mu) lepbar_v_pt.Rotate(smear_angle_mm_lepbar, transe_lepbar_axis);
	}
	
	lep_v_pt.Rotate(r.Rndm()*2.*TMath::Pi(), lep_v_pt_1);
	lepbar_v_pt.Rotate(r.Rndm()*2.*TMath::Pi(), lepbar_v_pt_1);
	
	lep.SetPxPyPzE(   lep_v_pt[0],    lep_v_pt[1],    lep_v_pt[2],    lep_pt_smear.E());
	lepbar.SetPxPyPzE(lepbar_v_pt[0], lepbar_v_pt[1], lepbar_v_pt[2], lepbar_pt_smear.E());

	if (applySmearing) {
	  smear_scale_jet_0 = h_pT_jet_smear->GetRandom();
	  smear_scale_jet_1 = h_pT_jet_smear->GetRandom();
	}
	
	if (i_jets==0) {
	  jetbar_pt_smear = tlv_jetbar;
	  jetbar_pt_smear.SetPtEtaPhiE(jetbar_pt_smear.Pt()*smear_scale_jet_0,jetbar_pt_smear.Eta(),jetbar_pt_smear.Phi(),jetbar_pt_smear.E());
	  jet_pt_smear    = tlv_jet;
	  jet_pt_smear.SetPtEtaPhiE(jet_pt_smear.Pt()*smear_scale_jet_1,jet_pt_smear.Eta(),jet_pt_smear.Phi(),jet_pt_smear.E());
	  jetbar_nosmear  = tlv_jetbar;
	  jet_nosmear     = tlv_jet;
	}
	else {
	  jetbar_pt_smear = tlv_jet;
	  jetbar_pt_smear.SetPtEtaPhiE(jetbar_pt_smear.Pt()*smear_scale_jet_0,jetbar_pt_smear.Eta(),jetbar_pt_smear.Phi(),jetbar_pt_smear.E());
	  jet_pt_smear    = tlv_jetbar;
	  jet_pt_smear.SetPtEtaPhiE(jet_pt_smear.Pt()*smear_scale_jet_1,jet_pt_smear.Eta(),jet_pt_smear.Phi(),jet_pt_smear.E());
	  jetbar_nosmear  = tlv_jet;
	  jet_nosmear     = tlv_jetbar;
	}

	// let's define axis.transverse to jet momentun (for theta smearing accounting)
	TVector3 jet_v_pt_1    = jet_pt_smear.Vect();
	TVector3 jetbar_v_pt_1 = jetbar_pt_smear.Vect();
	
	TVector3 jet_axis = jet_v_pt_1;
	TVector3 transe_jet_axis = beam_axis1.Cross(jet_axis);
	
	TVector3 jetbar_axis = jetbar_v_pt_1;
	TVector3 transe_jetbar_axis = beam_axis1.Cross(jetbar_axis);
	
	TVector3 jet_v_pt = jet_v_pt_1;
	TVector3 jetbar_v_pt = jetbar_v_pt_1;

	Double_t smear_angle_jet    = 0.0;
	Double_t smear_angle_jetbar = 0.0;
	
	if (applySmearing) {
	  smear_angle_jet    = h_theta_jet_smear->GetRandom();
	  smear_angle_jetbar = h_theta_jetbar_smear->GetRandom();
	}

	jet_v_pt.Rotate(smear_angle_jet, transe_jet_axis);
	jetbar_v_pt.Rotate(smear_angle_jetbar, transe_jetbar_axis);
	
	jet_v_pt.Rotate(r.Rndm()*2.*TMath::Pi(), jet_v_pt_1);
	jetbar_v_pt.Rotate(r.Rndm()*2.*TMath::Pi(), jetbar_v_pt_1);
	
	jet.SetPxPyPzE(   jet_v_pt[0],    jet_v_pt[1],    jet_v_pt[2],    jet_pt_smear.E());
	jetbar.SetPxPyPzE(jetbar_v_pt[0], jetbar_v_pt[1], jetbar_v_pt[2], jetbar_pt_smear.E());
	
	TLorentzVector mlbarb, mlbbar;
	mlbarb = lepbar + jet;
	mlbbar = lep + jetbar;
	
	TVector3 jet_v = jet.Vect();
	TVector3 lep_v = lep.Vect();
	TVector3 jetbar_v = jetbar.Vect();
	TVector3 lepbar_v = lepbar.Vect();

	Double_t Emiss_x_nosmear = Emiss_x;
	Double_t Emiss_y_nosmear = Emiss_y;
	Double_t Emiss_x_smear, Emiss_y_smear;

	Emiss_x_smear = Emiss_x_nosmear + (lep_nosmear.Px() - lep.Px()) + (lepbar_nosmear.Px() - lepbar.Px()) +
	  (jet_nosmear.Px() - jet.Px()) + (jetbar_nosmear.Px() - jetbar.Px());
	Emiss_y_smear = Emiss_y_nosmear  + (lep_nosmear.Py() - lep.Py()) + (lepbar_nosmear.Py() - lepbar.Py()) +
	  (jet_nosmear.Py() - jet.Py()) + (jetbar_nosmear.Py() - jetbar.Py());
	
	
	if (print) {
	  std::cout << "  Smearing iteration " << i_smear << ": " << std::endl;
	  std::cout << "   Lepton scale smearing: " << smear_scale_0 << ", " << smear_scale_1 << std::endl;
	  std::cout << "   Lepton angle smearing: " << smear_angle_ee_lep << ", "
	  		    << smear_angle_ee_lepbar << ", " <<  smear_angle_mm_lep << ", " << smear_angle_mm_lepbar << std::endl;
	  std::cout << "   Jet scale smearing: " <<  smear_scale_jet_0 << ", " << smear_scale_jet_1 << std::endl;
	  std::cout << "   Jet angle smearing: " <<  smear_angle_jet << ", " <<  smear_angle_jetbar << std::endl;;
	  std::cout << "   Px before smear (l, lbar, j, jbar, met): " << lep_nosmear.Px() << ",  " << lepbar_nosmear.Px()
	  		    << ",  " << jet_nosmear.Px() << ",  " << jetbar_nosmear.Px() << ", " << Emiss_x_nosmear  << std::endl;
	  std::cout << "   Px after  smear (l, lbar, j, jbar, met): " << lep.Px() << ",  " << lepbar.Px() << ",  "
	  		    << jet.Px() << ",  " << jetbar.Px() << ", " << Emiss_x_smear  << std::endl;
	}

	
	Double_t a1 = (jet.E() + lepbar.E())*(m_W*m_W - m_lepbar*m_lepbar - m_nu*m_nu)-
	  lepbar.E()*(m_top*m_top - m_b*m_b - m_lepbar*m_lepbar - m_nu*m_nu)+
	  2.*jet.E()*lepbar.E()*lepbar.E() - 2.*lepbar.E()*(jet_v.Dot(lepbar_v));
	Double_t a2 = 2.*(jet.E()*lepbar.Px() - lepbar.E()*jet.Px());
	Double_t a3 = 2.*(jet.E()*lepbar.Py() - lepbar.E()*jet.Py());
	Double_t a4 = 2.*(jet.E()*lepbar.Pz() - lepbar.E()*jet.Pz());
	
	
	Double_t b1 = (jetbar.E() + lep.E())*(m_Wbar*m_Wbar - m_lep*m_lep - m_nubar*m_nubar)-
	  lep.E()*(m_topbar*m_topbar - m_bbar*m_bbar - m_lep*m_lep - m_nubar*m_nubar)+
	  2.*jetbar.E()*lep.E()*lep.E() - 2.*lep.E()*(jetbar_v.Dot(lep_v));
	Double_t b2 = 2.*(jetbar.E()*lep.Px() - lep.E()*jetbar.Px());
	Double_t b3 = 2.*(jetbar.E()*lep.Py() - lep.E()*jetbar.Py());
	Double_t b4 = 2.*(jetbar.E()*lep.Pz() - lep.E()*jetbar.Pz());
	
	
	Double_t c22 = (m_W*m_W - m_lepbar*m_lepbar - m_nu*m_nu)*(m_W*m_W - m_lepbar*m_lepbar - m_nu*m_nu) -
	  4.*(lepbar.E()*lepbar.E()-lepbar.Pz()*lepbar.Pz())*(a1/a4)*(a1/a4) -
	  4.*(m_W*m_W - m_lepbar*m_lepbar - m_nu*m_nu)*lepbar.Pz()*a1/a4;
	Double_t c21 = 4.*(m_W*m_W - m_lepbar*m_lepbar - m_nu*m_nu)*(lepbar.Px() - lepbar.Pz()*a2/a4) -
	  8.*(lepbar.E()*lepbar.E() - lepbar.Pz()*lepbar.Pz())*a1*a2/(a4*a4) -
	  8*lepbar.Px()*lepbar.Pz()*a1/a4;
	Double_t c20 = -4.*(lepbar.E()*lepbar.E() - lepbar.Px()*lepbar.Px()) -
	  4.*(lepbar.E()*lepbar.E() - lepbar.Pz()*lepbar.Pz())*(a2/a4)*(a2/a4) -
	  8.*lepbar.Px()*lepbar.Pz()*a2/a4;
	
	Double_t c11 = 4.*(m_W*m_W - m_lepbar*m_lepbar - m_nu*m_nu)*(lepbar.Py() - lepbar.Pz()*a3/a4) -
	  8.*(lepbar.E()*lepbar.E() - lepbar.Pz()*lepbar.Pz())*a1*a3/(a4*a4) -
	  8.*lepbar.Py()*lepbar.Pz()*a1/a4;
	Double_t c10 = -8.*(lepbar.E()*lepbar.E() - lepbar.Pz()*lepbar.Pz())*a2*a3/(a4*a4) +
	  8.*lepbar.Px()*lepbar.Py() - 8.*lepbar.Px()*lepbar.Pz()*a3/a4 -
	  8.*lepbar.Py()*lepbar.Pz()*a2/a4;
	Double_t c00 = -4.*(lepbar.E()*lepbar.E() - lepbar.Py()*lepbar.Py()) -
	  4.*(lepbar.E()*lepbar.E() - lepbar.Pz()*lepbar.Pz())*(a3/a4)*(a3/a4) -
	  8.*lepbar.Py()*lepbar.Pz()*a3/a4;
	
	//new norm	
	c22 = c22*a4*a4;
	c21 = c21*a4*a4;
	c20 = c20*a4*a4;
	c11 = c11*a4*a4;
	c10 = c10*a4*a4;
	c00 = c00*a4*a4;
	
	Double_t dp22 = (m_Wbar*m_Wbar - m_lep*m_lep - m_nubar*m_nubar)*(m_Wbar*m_Wbar - m_lep*m_lep  - m_nubar*m_nubar) -
	  4.*(lep.E()*lep.E() - lep.Pz()*lep.Pz())*(b1/b4)*(b1/b4) -
	  4.*(m_Wbar*m_Wbar - m_lep*m_lep - m_nubar*m_nubar)*lep.Pz()*b1/b4;
	Double_t dp21 = 4.*(m_Wbar*m_Wbar - m_lep*m_lep - m_nubar*m_nubar)*(lep.Px() - lep.Pz()*b2/b4) -
	  8.*(lep.E()*lep.E() - lep.Pz()*lep.Pz())*b1*b2/(b4*b4) - 8.*lep.Px()*lep.Pz()*b1/b4;
	Double_t dp20 = -4.*(lep.E()*lep.E() - lep.Px()*lep.Px()) -
	  4.*(lep.E()*lep.E() - lep.Pz()*lep.Pz())*(b2/b4)*(b2/b4)-
	  8.*lep.Px()*lep.Pz()*b2/b4;
	
	Double_t dp11 = 4*(m_Wbar*m_Wbar - m_lep*m_lep - m_nubar*m_nubar)*(lep.Py() - lep.Pz()*b3/b4) -
	  8.*(lep.E()*lep.E() - lep.Pz()*lep.Pz())*b1*b3/(b4*b4) -
	  8.*lep.Py()*lep.Pz()*b1/b4;
	Double_t dp10 = -8.*(lep.E()*lep.E() - lep.Pz()*lep.Pz())*b2*b3/(b4*b4) +
	  8.*lep.Px()*lep.Py()- 8.*lep.Px()*lep.Pz()*b3/b4 - 8.*lep.Py()*lep.Pz()*b2/b4;
	Double_t dp00 = -4.*(lep.E()*lep.E() - lep.Py()*lep.Py()) -
	  4.*(lep.E()*lep.E() - lep.Pz()*lep.Pz())*(b3/b4)*(b3/b4)-
	  8.*lep.Py()*lep.Pz()*b3/b4;
	
	Double_t d22 = dp22 + Emiss_x_smear * Emiss_x_smear * dp20 + Emiss_y_smear * Emiss_y_smear * dp00 +
	  Emiss_x_smear * Emiss_y_smear * dp10 + Emiss_x_smear * dp21 + Emiss_y_smear * dp11;
	Double_t d21 = -1. * dp21 - 2.* Emiss_x_smear * dp20 - Emiss_y_smear * dp10;
	Double_t d20 = dp20;
	Double_t d11 = -1. * dp11 - 2. * Emiss_y_smear * dp00 - Emiss_x_smear * dp10;
	Double_t d10 = dp10;
	Double_t d00 = dp00;
	
	d22 = d22*b4*b4;
	d21 = d21*b4*b4;
	d20 = d20*b4*b4;
	d11 = d11*b4*b4;
	d10 = d10*b4*b4;
	d00 = d00*b4*b4;
	
	Double_t h4 = c00*c00*d22*d22 + c11*d22*(c11*d00 - c00*d11) +
	  c00*c22*(d11*d11 - 2.*d00*d22) +
	  c22*d00*(c22*d00 - c11*d11);
	
	Double_t h3 = c00*d21*(2.*c00*d22 - c11*d11) + c00*d11*(2.*c22*d10 + c21*d11)+
	  c22*d00*(2.*c21*d00 - c11*d10) - c00*d22*(c11*d10 + c10*d11) -
	  2.*c00*d00*(c22*d21 + c21*d22) - d00*d11*(c11*c21 + c10*c22) +
	  c11*d00*(c11*d21 + 2*c10*d22);
	
	Double_t h2 = c00*c00*(2.*d22*d20 + d21*d21) - c00*d21*(c11*d10 + c10*d11)+
	  c11*d20*(c11*d00 - c00*d11) + c00*d10*(c22*d10 - c10*d22) +
	  c00*d11*(2.*c21*d10 + c20*d11) + (2.*c22*c20 + c21*c21)*d00*d00 -
	  2.*c00*d00*(c22*d20 + c21*d21 + c20*d22)+
	  c10*d00*(2.*c11*d21 + c10*d22) - d00*d10*(c11*c21 + c10*c22)-
	  d00*d11*(c11*c20 + c10*c21);
	
	Double_t h1 = c00*d21*(2.*c00*d20 - c10*d10) - c00*d20*(c11*d10 + c10*d11) +
	  c00*d10*(c21*d10 + 2*c20*d11) - 2*c00*d00*(c21*d20 + c20*d21) +
	  c10*d00*(2.*c11*d20 + c10*d21) + c20*d00*(2*c21*d00 - c10*d11) -
	  d00*d10*(c11*c20 + c10*c21);
	
	Double_t h0 = c00*c00*d20*d20 + c10*d20*(c10*d00 - c00*d10) +
	  c20*d10*(c00*d10 - c10*d00) + c20*d00*(c20*d00 - 2.*c00*d20);
	
	ROOT::Math::Polynomial poly = ROOT::Math::Polynomial(h0,h1,h2,h3,h4);
	std::vector<Double_t> roots = poly.FindRealRoots();
	
	n_solutions = roots.size();
	if (print) {
	  std::cout << "   number of polynom roots:" << n_solutions << std::endl;
	}
	// ---> beginning of solution selection
	
	// Double_t mtt_min_diff = 100000.; unused!
	Double_t mtt_min = 100000.;
	
	Double_t Top_reco_px = 10000.;
	Double_t Top_reco_py = 10000.;
	Double_t Top_reco_pz = 10000.;
	
	Double_t Topbar_reco_px = 10000.;
	Double_t Topbar_reco_py = 10000.;
	Double_t Topbar_reco_pz = 10000.;

	if (n_solutions==0) continue;
	nIterations++;
	
	for (std::vector<Double_t>::iterator it = roots.begin() ; it!=roots.end() ; ++it){
	  p_nu_x_close =  *it;
	  Double_t c0 = c00;
	  Double_t c1 = c11 + c10*p_nu_x_close;
	  Double_t c2 = c22 + c21*p_nu_x_close + c20*p_nu_x_close*p_nu_x_close;
	  
	  Double_t d0 = d00;
	  Double_t d1 = d11 + d10*p_nu_x_close;
	  Double_t d2 = d22 + d21*p_nu_x_close + d20*p_nu_x_close*p_nu_x_close;
	  
	  p_nu_y = (c0*d2 -c2*d0)/(c1*d0 - c0*d1);
	  p_nu_z = (-1.*a1 -a2*p_nu_x_close - a3*p_nu_y)/a4;
	  delta_z = nu.Pz()-p_nu_z;
	  p_nubar_x = Emiss_x_smear - p_nu_x_close;
	  p_nubar_y = Emiss_y_smear - p_nu_y;
	  p_nubar_z = (-1.*b1 -b2*p_nubar_x- b3*p_nubar_y)/b4;
	  Double_t E_nu_calc = sqrt(p_nu_x_close*p_nu_x_close + p_nu_y*p_nu_y + p_nu_z*p_nu_z);
	  nu_calc.SetPxPyPzE(p_nu_x_close,p_nu_y,p_nu_z,E_nu_calc);
	  
	  Double_t E_nubar_calc = sqrt(p_nubar_x*p_nubar_x + p_nubar_y*p_nubar_y + p_nubar_z*p_nubar_z);
	  nubar_calc.SetPxPyPzE(p_nubar_x,p_nubar_y,p_nubar_z,E_nubar_calc);
	  
	  Top1 = nu_calc + lepbar + jet;
	  Top2 = nu + lepbar + jet;
	  
	  Top_calc    = nu_calc    + lepbar + jet;
	  Topbar_calc = nubar_calc + lep    + jetbar;
	  
	  TLorentzVector top_antitop_r, top_antitop_t;
	  top_antitop_t = top + topbar;
	  top_antitop_r = Top_calc + Topbar_calc;
	  
	  mtt_val.push_back(top_antitop_r.M());
	  
	  // mtt_truth = top_antitop_t.M(); unused!
	  
	  if(top_antitop_r.M() < mtt_min ) {
	    mtt_min = top_antitop_r.M();
	    
	    Top_mass_f = Top_calc.M();
	    Topbar_mass_f = Topbar_calc.M();
	    Top_mass_f_1 = Top2.M();
	    
	    Top_reco_px = Top_calc.Px();
	    Top_reco_py = Top_calc.Py();
	    Top_reco_pz = Top_calc.Pz();
	    
	    Topbar_reco_px = Topbar_calc.Px();
	    Topbar_reco_py = Topbar_calc.Py();
	    Topbar_reco_pz = Topbar_calc.Pz();
	    
	    p_nu_x_close_f = *it;
	    p_nu_y_f = p_nu_y;
	    p_nu_z_f = p_nu_z;
	    
	    p_nubar_x_f = p_nubar_x;
	    p_nubar_y_f = p_nubar_y;
	    p_nubar_z_f = p_nubar_z;
	    
	    delta_TOP_mass = Top_calc.M() - m_top;
	  } // ---> end of solution selection
	} // ---> end of the loop over solutions      
	
	Double_t weight_s_i1 =1.;
	Double_t weight_s_i2 =1.;
	
	TVector3 top_p_i(Top_reco_px,Top_reco_py,Top_reco_pz);
	TVector3 topbar_p_i(Topbar_reco_px,Topbar_reco_py,Topbar_reco_pz);
	
	TAxis *xaxis1 = h_mlblb_smear->GetXaxis();
	Int_t binx1 = xaxis1->FindBin(mlbarb.M());
	weight_s_i1 = h_mlblb_smear->GetBinContent(binx1);
	
	TAxis *xaxis2 = h_mlblb_smear->GetXaxis();
	Int_t binx2 = xaxis2->FindBin(mlbbar.M());
	weight_s_i2 = h_mlblb_smear->GetBinContent(binx2);
	
	Vec_Top[i_jets][0] = Vec_Top[i_jets][0] + weight_s_i1*weight_s_i2*Top_reco_px;
	Vec_Top[i_jets][1] = Vec_Top[i_jets][1] + weight_s_i1*weight_s_i2*Top_reco_py;
	Vec_Top[i_jets][2] = Vec_Top[i_jets][2] + weight_s_i1*weight_s_i2*Top_reco_pz;
		
	Vec_Topbar[i_jets][0] = Vec_Topbar[i_jets][0] + weight_s_i1*weight_s_i2*Topbar_reco_px;
	Vec_Topbar[i_jets][1] = Vec_Topbar[i_jets][1] + weight_s_i1*weight_s_i2*Topbar_reco_py;
	Vec_Topbar[i_jets][2] = Vec_Topbar[i_jets][2] + weight_s_i1*weight_s_i2*Topbar_reco_pz;
	
	weight_s_sum = weight_s_sum + weight_s_i1*weight_s_i2;
	
	i_analyzed_event = 1;

	if (print) {
	  std::cout << "   weight            : " << weight_s_sum << std::endl;
	  std::cout << "   (px, py, pz)[t]   : " << Top_reco_px << ", " << Top_reco_py << ", " << Top_reco_pz << std::endl;
	  std::cout << "   (px, py, pz)[tbar]: " << Topbar_reco_px << ", " << Topbar_reco_py << ", " << Topbar_reco_pz << std::endl;
	}	
	
      } // end loop over smear variations

      if (print)
	std::cout << " Number of iteration with solutions : " << nIterations << std::endl;
      
      weights_com.push_back(weight_s_sum);
      
      if(i_jets == 0) {lep0 = lep; lepbar0 = lepbar; }
      if(i_jets == 1) {lep1 = lep; lepbar1 = lepbar; }
      
    } // end of loop over jets                                 

    // if not usefil solutions -> go to the next event
    if(i_analyzed_event == 0 || weights_com.size() == 0) {
      return return_top_tbar;
    }
    
    
    if (print)
      std::cout << " Weight sum of jet combinatorics (0, 1): " << weights_com[0] << ", " << weights_com[1] << std::endl;
      
   
    TVector3 top_p_sum;
    TVector3 topbar_p_sum;
    Int_t i_weight_sel = 0;
    
    if(weights_com[0] > weights_com[1] ) {
      if(weights_com[0] > 0. && Vec_Top[0][0] !=0) {

	top_p_sum.SetX(Vec_Top[0][0]/weights_com[0]);
	top_p_sum.SetY(Vec_Top[0][1]/weights_com[0]);
	top_p_sum.SetZ(Vec_Top[0][2]/weights_com[0]);
	
	topbar_p_sum.SetX(Vec_Topbar[0][0]/weights_com[0]);
	topbar_p_sum.SetY(Vec_Topbar[0][1]/weights_com[0]);
	topbar_p_sum.SetZ(Vec_Topbar[0][2]/weights_com[0]);
	
	lplus  = lepbar0;
	lminus = lep0;
	
	i_weight_sel = 1;
	
      }
      
    } else {
      if(weights_com[1] > 0. && Vec_Top[1][0] !=0) {
	
	top_p_sum.SetX(Vec_Top[1][0]/weights_com[1]);
	top_p_sum.SetY(Vec_Top[1][1]/weights_com[1]);
	top_p_sum.SetZ(Vec_Top[1][2]/weights_com[1]);
	
	topbar_p_sum.SetX(Vec_Topbar[1][0]/weights_com[1]);
	topbar_p_sum.SetY(Vec_Topbar[1][1]/weights_com[1]);
	topbar_p_sum.SetZ(Vec_Topbar[1][2]/weights_com[1]);
	
	lplus  = lepbar1;
	lminus = lep1;
	
	i_weight_sel = 1;
	
      }
    }
    
    if(i_weight_sel == 0) {
      return return_top_tbar;
    }
    
    TLorentzVector Top_fin;
    TLorentzVector Topbar_fin;
    
    Double_t Top_fin_M =  m_top;
    Double_t Topbar_fin_M = m_topbar; 
    
    Double_t Top_fin_E    = sqrt(top_p_sum.Mag2()    + Top_fin_M*Top_fin_M);
    Double_t Topbar_fin_E = sqrt(topbar_p_sum.Mag2() + Topbar_fin_M*Topbar_fin_M);
    
    Top_fin.SetPxPyPzE(top_p_sum[0],    top_p_sum[1],    top_p_sum[2],    Top_fin_E);
    Topbar_fin.SetPxPyPzE(topbar_p_sum[0], topbar_p_sum[1], topbar_p_sum[2], Topbar_fin_E);
    
    return_top_tbar[0] = Top_fin;
    return_top_tbar[1] = Topbar_fin;
    return return_top_tbar;
  } 
};
