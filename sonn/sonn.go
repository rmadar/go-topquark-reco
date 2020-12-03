package sonn

import (
	"fmt"
	"log"
	"math"

	"go-hep.org/x/exp/roots"
	"go-hep.org/x/hep/fmom"
	"go-hep.org/x/hep/groot"
	"go-hep.org/x/hep/groot/rhist"
	"go-hep.org/x/hep/hbook"
	"go-hep.org/x/hep/hbook/rootcnv"
	"golang.org/x/exp/rand"
	"gonum.org/v1/gonum/spatial/r3"
)

// SmearingHistos holds the histograms used to smear 4-vectors.
type SmearingHistos struct {
	PtLepEE    hbook.Rand1D
	PtLepMu    hbook.Rand1D
	ThetaLepEE hbook.Rand1D
	ThetaLepMu hbook.Rand1D

	PtJet       hbook.Rand1D
	ThetaJet    hbook.Rand1D
	ThetaJetBar hbook.Rand1D

	Mlblb *hbook.H1D
}

// NewSmearingHistos returns a set of smearing histograms from the provided
// ROOT file name and a seed for the histogram distributions.
func NewSmearingHistos(fname string, seed uint64) (*SmearingHistos, error) {
	f, err := groot.Open(fname)
	if err != nil {
		return nil, fmt.Errorf("could not open ROOT file %q: %w", fname, err)
	}
	defer f.Close()

	get := func(name string) *hbook.H1D {
		o, e := f.Get(name)
		if e != nil && err != nil {
			err = fmt.Errorf("could not load histo %q: %w", name, e)
			return nil
		}
		return rootcnv.H1D(o.(*rhist.H1F))
	}

	gethr := func(name string) hbook.Rand1D {
		h := get(name)
		if h == nil {
			return hbook.Rand1D{}
		}
		return hbook.NewRand1D(h, rand.NewSource(seed))
	}

	sh := &SmearingHistos{
		PtLepEE:     gethr("h_pT_lep_ee_smear"),
		PtLepMu:     gethr("h_pT_lep_mm_smear"),
		Mlblb:       get("mlblb_smear"),
		ThetaJet:    gethr("h_theta_jet_smear"),
		ThetaJetBar: gethr("h_theta_jetbar_smear"),
		PtJet:       gethr("h_pT_lep_mm_smear"),
		ThetaLepEE:  gethr("h_theta_lep_ee_smear"),
		ThetaLepMu:  gethr("h_theta_lep_mm_smear"),
	}

	return sh, err
}

func Sonnenschein(
	lepTLV, lepbarTLV fmom.P4, pdgIDLep, pdgIDLepBar int,
	jetTLV, jetbarTLV fmom.P4, isbJet, isbJetbar bool,
	emissx, emissy float64,
	smearHs *SmearingHistos,
) []fmom.PxPyPzE {

	// debug.
	const debug = true

	// Smearing parameters
	var Nsmear = 10
	const applySmearing = true
	if !applySmearing {
		Nsmear = 1
	}

	// run the reconstruction only if the two jets are btags.
	if !isbJet || !isbJetbar {
		return nil
	}

	var (
		lepIsE     = iabs(pdgIDLep) == 11
		lepIsMu    = iabs(pdgIDLep) == 13
		lepBarIsE  = iabs(pdgIDLepBar) == 11
		lepBarIsMu = iabs(pdgIDLepBar) == 13
	)

	// particle masses in GeV
	const (
		m_W      = 80.379
		m_Wbar   = 80.379
		m_b      = 4.18
		m_bbar   = 4.18
		m_nu     = 0.0
		m_nubar  = 0.0
		m_top    = 172.5
		m_topbar = 172.5
		// m_e      = 0.000511; unused!
		// m_mu     = 0.10566; unused!
		// m_tau    = 1.77686; unused!
		m_lep    = 0.0
		m_lepbar = 0.0
	)

	var (
		nu                            fmom.PxPyPzE
		jet, jetbar                   fmom.PxPyPzE
		lep, lepbar                   fmom.PxPyPzE
		lep_pt_smear, lepbar_pt_smear fmom.PxPyPzE
		lep_nosmear, lepbar_nosmear   fmom.PxPyPzE
		jet_pt_smear, jetbar_pt_smear fmom.PxPyPzE
		jet_nosmear, jetbar_nosmear   fmom.PxPyPzE
		top, topbar                   fmom.PxPyPzE
		Top_calc, Topbar_calc         fmom.PxPyPzE
		nu_calc, nubar_calc           fmom.PxPyPzE
		Top1, Top2                    fmom.PxPyPzE

		Vec_Top    [2][3]float64
		Vec_Topbar [2][3]float64

		weights_com []float64

		i_analyzed_event int
		n_solutions      int

		p_nu_y = -999.0
		p_nu_z = -999.0

		p_nubar_x = -999.0
		p_nubar_y = -999.0
		p_nubar_z = -999.0

		// delta_z = 100000.0

		p_nu_x_close = 100000.0

		// p_nu_x_close_f = -999.0
		//	p_nu_y_f       = -999.0
		//	p_nu_z_f       = -999.0
		//	p_nubar_x_f    = -999.0
		//	p_nubar_y_f    = -999.0
		//	p_nubar_z_f    = -999.0
		//	delta_TOP_mass = 100000.0
		//	Top_mass_f     = 100000.0
		//	Topbar_mass_f  = 100000.0
		//	Top_mass_f_1   = 100000.0

		mtt_val []float64
	)

	// loop over jets
	for i_jets := 0; i_jets < 2; i_jets++ {
		var (
			weight_s_sum float64
			rnd          = rand.New(rand.NewSource(1))
		)

		if debug {
			log.Printf(" Jets combination %d:", i_jets)
		}

		nIterations := 0
		for i_smear := 0; i_smear < Nsmear; i_smear++ {
			var (
				smear_scale_0     = 1.0
				smear_scale_1     = 1.0
				smear_scale_jet_0 = 1.0
				smear_scale_jet_1 = 1.0
			)

			if applySmearing {
				switch {
				case lepIsE:
					smear_scale_0 = smearHs.PtLepEE.Rand()
				case lepIsMu:
					smear_scale_0 = smearHs.PtLepMu.Rand()
				default:
					smear_scale_0 = 1
				}

				switch {
				case lepBarIsE:
					smear_scale_1 = smearHs.PtLepEE.Rand()
				case lepBarIsMu:
					smear_scale_1 = smearHs.PtLepMu.Rand()
				default:
					smear_scale_1 = 1
				}
			}

			lepbar_pt_smear = newPtEtaPhiM(
				lepbarTLV.Pt()*smear_scale_0,
				lepbarTLV.Eta(),
				lepbarTLV.Phi(),
				m_lepbar,
			)

			lep_pt_smear = newPtEtaPhiM(
				lepTLV.Pt()*smear_scale_1,
				lepTLV.Eta(),
				lepTLV.Phi(),
				m_lepbar,
			)

			lepbar_nosmear.Set(lepbarTLV)
			lep_nosmear.Set(lepTLV)

			// let's define the axis transverse to lepton momentum
			// for theta smearing accounting.
			var (
				lep_v_pt_1 = r3.Vec{
					X: lep_pt_smear.Px(),
					Y: lep_pt_smear.Py(),
					Z: lep_pt_smear.Pz(),
				}
				lepbar_v_pt_1 = r3.Vec{
					X: lepbar_pt_smear.Px(),
					Y: lepbar_pt_smear.Py(),
					Z: lepbar_pt_smear.Pz(),
				}
				beam_axis1 = r3.Vec{Z: 1}

				lep_axis        = lep_v_pt_1
				transe_lep_axis = beam_axis1.Cross(lep_axis)

				lepbar_axis        = lepbar_v_pt_1
				transe_lepbar_axis = beam_axis1.Cross(lepbar_axis)

				lep_v_pt    = lep_v_pt_1
				lepbar_v_pt = lepbar_v_pt_1

				smear_angle_ee_lep    = 0.0
				smear_angle_ee_lepbar = 0.0
				smear_angle_mm_lep    = 0.0
				smear_angle_mm_lepbar = 0.0
			)

			if applySmearing {
				smear_angle_ee_lep = smearHs.ThetaLepEE.Rand()
				smear_angle_ee_lepbar = smearHs.ThetaLepEE.Rand()
				smear_angle_mm_lep = smearHs.ThetaLepMu.Rand()
				smear_angle_mm_lepbar = smearHs.ThetaLepMu.Rand()
			}

			switch {
			case lepTLV.Pt() > lepbarTLV.Pt():
				switch {
				case lepIsE:
					lep_v_pt = lep_v_pt.Rotate(
						smear_angle_ee_lep,
						transe_lep_axis,
					)
				case lepIsMu:
					lep_v_pt = lep_v_pt.Rotate(
						smear_angle_mm_lep,
						transe_lep_axis,
					)
				}
			default:
				switch {
				case lepBarIsE:
					lepbar_v_pt = lepbar_v_pt.Rotate(
						smear_angle_ee_lepbar,
						transe_lepbar_axis,
					)
				case lepBarIsMu:
					lepbar_v_pt = lepbar_v_pt.Rotate(
						smear_angle_mm_lepbar,
						transe_lepbar_axis,
					)
				}
			}

			lep_v_pt = lep_v_pt.Rotate(
				rnd.Float64()*2*math.Pi,
				lep_v_pt_1,
			)
			lepbar_v_pt = lepbar_v_pt.Rotate(
				rnd.Float64()*2*math.Pi,
				lepbar_v_pt_1,
			)

			lep = fmom.NewPxPyPzE(
				lep_v_pt.X, lep_v_pt.Y, lep_v_pt.Z,
				lep_pt_smear.E(),
			)
			lepbar = fmom.NewPxPyPzE(
				lepbar_v_pt.X, lepbar_v_pt.Y, lepbar_v_pt.Z,
				lepbar_pt_smear.E(),
			)

			if applySmearing {
				smear_scale_jet_0 = smearHs.PtJet.Rand()
				smear_scale_jet_1 = smearHs.PtJet.Rand()
			}

			switch i_jets {
			case 0:
				jetbar_pt_smear = newPtEtaPhiE(
					jetbarTLV.Pt()*smear_scale_jet_0,
					jetbarTLV.Eta(),
					jetbarTLV.Phi(),
					jetbarTLV.E(),
				)

				jet_pt_smear = newPtEtaPhiE(
					jetTLV.Pt()*smear_scale_jet_1,
					jetTLV.Eta(),
					jetTLV.Phi(),
					jetTLV.E(),
				)

				jetbar_nosmear.Set(jetbarTLV)
				jet_nosmear.Set(jetTLV)

			case 1:
				jetbar_pt_smear = newPtEtaPhiE(
					jetTLV.Pt()*smear_scale_jet_0,
					jetTLV.Eta(),
					jetTLV.Phi(),
					jetTLV.E(),
				)

				jet_pt_smear = newPtEtaPhiE(
					jetbarTLV.Pt()*smear_scale_jet_1,
					jetbarTLV.Eta(),
					jetbarTLV.Phi(),
					jetbarTLV.E(),
				)

				jetbar_nosmear.Set(jetTLV)
				jet_nosmear.Set(jetbarTLV)
			}

			// let's define the transverse axis to the jet momentum
			// for theta smearing accounting.
			var (
				jet_v_pt_1 = r3.Vec{
					X: jet_pt_smear.Px(),
					Y: jet_pt_smear.Py(),
					Z: jet_pt_smear.Pz(),
				}
				jetbar_v_pt_1 = r3.Vec{
					X: jetbar_pt_smear.Px(),
					Y: jetbar_pt_smear.Py(),
					Z: jetbar_pt_smear.Pz(),
				}

				jet_axis        = jet_v_pt_1
				transe_jet_axis = beam_axis1.Cross(jet_axis)

				jetbar_axis        = jetbar_v_pt_1
				transe_jetbar_axis = beam_axis1.Cross(jetbar_axis)

				jet_v_pt    = jet_v_pt_1
				jetbar_v_pt = jetbar_v_pt_1

				smear_angle_jet    = 0.0
				smear_angle_jetbar = 0.0
			)

			if applySmearing {
				smear_angle_jet = smearHs.ThetaJet.Rand()
				smear_angle_jetbar = smearHs.ThetaJetBar.Rand()
			}

			jet_v_pt = jet_v_pt.Rotate(smear_angle_jet, transe_jet_axis)
			jetbar_v_pt = jetbar_v_pt.Rotate(smear_angle_jetbar, transe_jetbar_axis)

			jet_v_pt = jet_v_pt.Rotate(rnd.Float64()*2*math.Pi, jet_v_pt_1)
			jetbar_v_pt = jetbar_v_pt.Rotate(rnd.Float64()*2*math.Pi, jetbar_v_pt_1)

			jet = fmom.NewPxPyPzE(
				jet_v_pt.X,
				jet_v_pt.Y,
				jet_v_pt.Z,
				jet_pt_smear.E(),
			)
			jetbar = fmom.NewPxPyPzE(
				jetbar_v_pt.X,
				jetbar_v_pt.Y,
				jetbar_v_pt.Z,
				jetbar_pt_smear.E(),
			)

			var (
				mlbarb fmom.PxPyPzE
				mlbbar fmom.PxPyPzE
			)
			mlbarb.Set(fmom.Add(&lepbar, &jet))
			mlbbar.Set(fmom.Add(&lep, &jetbar))

			var (
				jet_v    = r3.Vec{X: jet.Px(), Y: jet.Py(), Z: jet.Pz()}
				lep_v    = r3.Vec{X: lep.Px(), Y: lep.Py(), Z: lep.Pz()}
				jetbar_v = r3.Vec{X: jetbar.Px(), Y: jetbar.Py(), Z: jetbar.Pz()}
				lepbar_v = r3.Vec{X: lepbar.Px(), Y: lepbar.Py(), Z: lepbar.Pz()}

				Emiss_x_nosmear = emissx
				Emiss_y_nosmear = emissy
				Emiss_x_smear   float64
				Emiss_y_smear   float64
			)

			Emiss_x_smear = Emiss_x_nosmear +
				(lep_nosmear.Px() - lep.Px()) +
				(lepbar_nosmear.Px() - lepbar.Px()) +
				(jet_nosmear.Px() - jet.Px()) +
				(jetbar_nosmear.Px() - jetbar.Px())
			Emiss_y_smear = Emiss_y_nosmear +
				(lep_nosmear.Py() - lep.Py()) +
				(lepbar_nosmear.Py() - lepbar.Py()) +
				(jet_nosmear.Py() - jet.Py()) +
				(jetbar_nosmear.Py() - jetbar.Py())

			if debug {
				log.Printf("  Smearing iteration %d:", i_smear)
				log.Printf("   Lepton scale smearing: %g, %g", smear_scale_0, smear_scale_1)
				log.Printf("   Lepton angle smearing: %g, %g, %g, %g",
					smear_angle_ee_lep,
					smear_angle_ee_lepbar,
					smear_angle_mm_lep,
					smear_angle_mm_lepbar,
				)
				log.Printf("   Jet scale smearing: %g, %g",
					smear_scale_jet_0,
					smear_scale_jet_1,
				)
				log.Printf("   Jet angle smearing: %g, %g",
					smear_angle_jet,
					smear_angle_jetbar,
				)
				log.Printf("   Px before smear (l, lbar, j, jbar, met): %g, %g, %g, %g, %g",
					lep_nosmear.Px(),
					lepbar_nosmear.Px(),
					jet_nosmear.Px(),
					jetbar_nosmear.Px(),
					Emiss_x_nosmear,
				)
				log.Printf("   Px after  smear (l, lbar, j, jbar, met): %g, %g, %g, %g, %g",
					lep.Px(),
					lepbar.Px(),
					jet.Px(),
					jetbar.Px(),
					Emiss_x_smear,
				)
			}

			a1 := (jet.E()+lepbar.E())*(m_W*m_W-m_lepbar*m_lepbar-m_nu*m_nu) -
				lepbar.E()*(m_top*m_top-m_b*m_b-m_lepbar*m_lepbar-m_nu*m_nu) +
				2*jet.E()*lepbar.E()*lepbar.E() - 2*lepbar.E()*(jet_v.Dot(lepbar_v))
			a2 := 2 * (jet.E()*lepbar.Px() - lepbar.E()*jet.Px())
			a3 := 2 * (jet.E()*lepbar.Py() - lepbar.E()*jet.Py())
			a4 := 2 * (jet.E()*lepbar.Pz() - lepbar.E()*jet.Pz())

			b1 := (jetbar.E()+lep.E())*(m_Wbar*m_Wbar-m_lep*m_lep-m_nubar*m_nubar) -
				lep.E()*(m_topbar*m_topbar-m_bbar*m_bbar-m_lep*m_lep-m_nubar*m_nubar) +
				2.*jetbar.E()*lep.E()*lep.E() - 2.*lep.E()*(jetbar_v.Dot(lep_v))
			b2 := 2. * (jetbar.E()*lep.Px() - lep.E()*jetbar.Px())
			b3 := 2. * (jetbar.E()*lep.Py() - lep.E()*jetbar.Py())
			b4 := 2. * (jetbar.E()*lep.Pz() - lep.E()*jetbar.Pz())

			c22 := (m_W*m_W-m_lepbar*m_lepbar-m_nu*m_nu)*(m_W*m_W-m_lepbar*m_lepbar-m_nu*m_nu) -
				4.*(lepbar.E()*lepbar.E()-lepbar.Pz()*lepbar.Pz())*(a1/a4)*(a1/a4) -
				4.*(m_W*m_W-m_lepbar*m_lepbar-m_nu*m_nu)*lepbar.Pz()*a1/a4
			c21 := 4.*(m_W*m_W-m_lepbar*m_lepbar-m_nu*m_nu)*(lepbar.Px()-lepbar.Pz()*a2/a4) -
				8.*(lepbar.E()*lepbar.E()-lepbar.Pz()*lepbar.Pz())*a1*a2/(a4*a4) -
				8*lepbar.Px()*lepbar.Pz()*a1/a4
			c20 := -4.*(lepbar.E()*lepbar.E()-lepbar.Px()*lepbar.Px()) -
				4.*(lepbar.E()*lepbar.E()-lepbar.Pz()*lepbar.Pz())*(a2/a4)*(a2/a4) -
				8.*lepbar.Px()*lepbar.Pz()*a2/a4

			c11 := 4.*(m_W*m_W-m_lepbar*m_lepbar-m_nu*m_nu)*(lepbar.Py()-lepbar.Pz()*a3/a4) -
				8.*(lepbar.E()*lepbar.E()-lepbar.Pz()*lepbar.Pz())*a1*a3/(a4*a4) -
				8.*lepbar.Py()*lepbar.Pz()*a1/a4
			c10 := -8.*(lepbar.E()*lepbar.E()-lepbar.Pz()*lepbar.Pz())*a2*a3/(a4*a4) +
				8.*lepbar.Px()*lepbar.Py() - 8.*lepbar.Px()*lepbar.Pz()*a3/a4 -
				8.*lepbar.Py()*lepbar.Pz()*a2/a4
			c00 := -4.*(lepbar.E()*lepbar.E()-lepbar.Py()*lepbar.Py()) -
				4.*(lepbar.E()*lepbar.E()-lepbar.Pz()*lepbar.Pz())*(a3/a4)*(a3/a4) -
				8.*lepbar.Py()*lepbar.Pz()*a3/a4

			// new norm
			n44 := a4 * a4
			c22 *= n44
			c21 *= n44
			c20 *= n44
			c11 *= n44
			c10 *= n44
			c00 *= n44

			dp22 := (m_Wbar*m_Wbar-m_lep*m_lep-m_nubar*m_nubar)*(m_Wbar*m_Wbar-m_lep*m_lep-m_nubar*m_nubar) -
				4.*(lep.E()*lep.E()-lep.Pz()*lep.Pz())*(b1/b4)*(b1/b4) -
				4.*(m_Wbar*m_Wbar-m_lep*m_lep-m_nubar*m_nubar)*lep.Pz()*b1/b4
			dp21 := 4.*(m_Wbar*m_Wbar-m_lep*m_lep-m_nubar*m_nubar)*(lep.Px()-lep.Pz()*b2/b4) -
				8.*(lep.E()*lep.E()-lep.Pz()*lep.Pz())*b1*b2/(b4*b4) - 8.*lep.Px()*lep.Pz()*b1/b4
			dp20 := -4.*(lep.E()*lep.E()-lep.Px()*lep.Px()) -
				4.*(lep.E()*lep.E()-lep.Pz()*lep.Pz())*(b2/b4)*(b2/b4) -
				8.*lep.Px()*lep.Pz()*b2/b4

			dp11 := 4*(m_Wbar*m_Wbar-m_lep*m_lep-m_nubar*m_nubar)*(lep.Py()-lep.Pz()*b3/b4) -
				8.*(lep.E()*lep.E()-lep.Pz()*lep.Pz())*b1*b3/(b4*b4) -
				8.*lep.Py()*lep.Pz()*b1/b4
			dp10 := -8.*(lep.E()*lep.E()-lep.Pz()*lep.Pz())*b2*b3/(b4*b4) +
				8.*lep.Px()*lep.Py() - 8.*lep.Px()*lep.Pz()*b3/b4 - 8.*lep.Py()*lep.Pz()*b2/b4
			dp00 := -4.*(lep.E()*lep.E()-lep.Py()*lep.Py()) -
				4.*(lep.E()*lep.E()-lep.Pz()*lep.Pz())*(b3/b4)*(b3/b4) -
				8.*lep.Py()*lep.Pz()*b3/b4

			d22 := dp22 +
				Emiss_x_smear*Emiss_x_smear*dp20 +
				Emiss_y_smear*Emiss_y_smear*dp00 +
				Emiss_x_smear*Emiss_y_smear*dp10 +
				Emiss_x_smear*dp21 +
				Emiss_y_smear*dp11
			d21 := -1.*dp21 - 2.*Emiss_x_smear*dp20 - Emiss_y_smear*dp10
			d20 := dp20
			d11 := -1.*dp11 - 2.*Emiss_y_smear*dp00 - Emiss_x_smear*dp10
			d10 := dp10
			d00 := dp00

			d22 = d22 * b4 * b4
			d21 = d21 * b4 * b4
			d20 = d20 * b4 * b4
			d11 = d11 * b4 * b4
			d10 = d10 * b4 * b4
			d00 = d00 * b4 * b4

			h4 := c00*c00*d22*d22 + c11*d22*(c11*d00-c00*d11) +
				c00*c22*(d11*d11-2.*d00*d22) +
				c22*d00*(c22*d00-c11*d11)

			h3 := c00*d21*(2.*c00*d22-c11*d11) + c00*d11*(2.*c22*d10+c21*d11) +
				c22*d00*(2.*c21*d00-c11*d10) - c00*d22*(c11*d10+c10*d11) -
				2.*c00*d00*(c22*d21+c21*d22) - d00*d11*(c11*c21+c10*c22) +
				c11*d00*(c11*d21+2*c10*d22)

			h2 := c00*c00*(2.*d22*d20+d21*d21) - c00*d21*(c11*d10+c10*d11) +
				c11*d20*(c11*d00-c00*d11) + c00*d10*(c22*d10-c10*d22) +
				c00*d11*(2.*c21*d10+c20*d11) + (2.*c22*c20+c21*c21)*d00*d00 -
				2.*c00*d00*(c22*d20+c21*d21+c20*d22) +
				c10*d00*(2.*c11*d21+c10*d22) - d00*d10*(c11*c21+c10*c22) -
				d00*d11*(c11*c20+c10*c21)

			h1 := c00*d21*(2.*c00*d20-c10*d10) - c00*d20*(c11*d10+c10*d11) +
				c00*d10*(c21*d10+2*c20*d11) - 2*c00*d00*(c21*d20+c20*d21) +
				c10*d00*(2.*c11*d20+c10*d21) + c20*d00*(2*c21*d00-c10*d11) -
				d00*d10*(c11*c20+c10*c21)

			h0 := c00*c00*d20*d20 + c10*d20*(c10*d00-c00*d10) +
				c20*d10*(c00*d10-c10*d00) + c20*d00*(c20*d00-2.*c00*d20)

			var zs [4]complex128
			zs[0], zs[1], zs[2], zs[3] = roots.Poly4(h0, h1, h2, h3, h4)

			roots := make([]float64, 0, 4)
			const ε = 1e-15
			for _, z := range zs {
				if math.Abs(imag(z)) < ε {
					roots = append(roots, real(z))
				}
			}
			//sort.Float64s(roots)
			n_solutions = len(roots)
			if debug {
				log.Printf("   number of polynom roots: %d", n_solutions)
			}

			if n_solutions == 0 {
				continue
			}
			nIterations++

			// select roots
			var (
				mtt_min = 100000.0

				Top_reco_px = 10000.0
				Top_reco_py = 10000.0
				Top_reco_pz = 10000.0

				Topbar_reco_px = 10000.0
				Topbar_reco_py = 10000.0
				Topbar_reco_pz = 10000.0
			)

			for _, v := range roots {
				p_nu_x_close = v

				c0 := c00
				c1 := c11 + c10*p_nu_x_close
				c2 := c22 + c21*p_nu_x_close + c20*p_nu_x_close*p_nu_x_close

				d0 := d00
				d1 := d11 + d10*p_nu_x_close
				d2 := d22 + d21*p_nu_x_close + d20*p_nu_x_close*p_nu_x_close

				p_nu_y = (c0*d2 - c2*d0) / (c1*d0 - c0*d1)
				p_nu_z = (-1.*a1 - a2*p_nu_x_close - a3*p_nu_y) / a4
				// delta_z = nu.Pz() - p_nu_z
				p_nubar_x = Emiss_x_smear - p_nu_x_close
				p_nubar_y = Emiss_y_smear - p_nu_y
				p_nubar_z = (-1.*b1 - b2*p_nubar_x - b3*p_nubar_y) / b4

				E_nu_calc := math.Sqrt(p_nu_x_close*p_nu_x_close + p_nu_y*p_nu_y + p_nu_z*p_nu_z)
				nu_calc = fmom.NewPxPyPzE(p_nu_x_close, p_nu_y, p_nu_z, E_nu_calc)

				E_nubar_calc := math.Sqrt(p_nubar_x*p_nubar_x + p_nubar_y*p_nubar_y + p_nubar_z*p_nubar_z)
				nubar_calc = fmom.NewPxPyPzE(p_nubar_x, p_nubar_y, p_nubar_z, E_nubar_calc)

				Top1.Set(fmom.Add(fmom.Add(&nu_calc, &lepbar), &jet))
				Top2.Set(fmom.Add(fmom.Add(&nu, &lepbar), &jet))

				Top_calc.Set(fmom.Add(fmom.Add(&nu_calc, &lepbar), &jet))
				Topbar_calc.Set(fmom.Add(fmom.Add(&nubar_calc, &lep), &jetbar))

				var (
					top_antitop_r fmom.PxPyPzE
					top_antitop_t fmom.PxPyPzE
				)
				top_antitop_t.Set(fmom.Add(&top, &topbar))
				top_antitop_r.Set(fmom.Add(&Top_calc, &Topbar_calc))

				mtt_val = append(mtt_val, top_antitop_r.M())

				if top_antitop_r.M() < mtt_min {
					mtt_min = top_antitop_r.M()

					//Top_mass_f = Top_calc.M()
					//Topbar_mass_f = Topbar_calc.M()
					//Top_mass_f_1 = Top2.M()

					Top_reco_px = Top_calc.Px()
					Top_reco_py = Top_calc.Py()
					Top_reco_pz = Top_calc.Pz()

					Topbar_reco_px = Topbar_calc.Px()
					Topbar_reco_py = Topbar_calc.Py()
					Topbar_reco_pz = Topbar_calc.Pz()

					//p_nu_x_close_f = v
					//p_nu_y_f = p_nu_y
					//p_nu_z_f = p_nu_z

					//p_nubar_x_f = p_nubar_x
					//p_nubar_y_f = p_nubar_y
					//p_nubar_z_f = p_nubar_z

					//delta_TOP_mass = Top_calc.M() - m_top
				} // ---> end of solution selection
			}

			var (
				// FIXME?
				// top_p_i    = r3.Vec{X: Top_reco_px, Y: Top_reco_py, Z: Top_reco_pz}
				// topbar_p_i = r3.Vec{X: Topbar_reco_px, Y: Topbar_reco_py, Z: Topbar_reco_pz}

				binx1 = hbook.Bin1Ds(smearHs.Mlblb.Binning.Bins).IndexOf(mlbarb.M())

				binx2       = hbook.Bin1Ds(smearHs.Mlblb.Binning.Bins).IndexOf(mlbbar.M())
				weight_s_i1 float64
				weight_s_i2 float64
			)

			switch binx1 {
			case hbook.UnderflowBin1D:
				weight_s_i1 = smearHs.Mlblb.Binning.Underflow().SumW()
			case hbook.OverflowBin1D:
				weight_s_i1 = smearHs.Mlblb.Binning.Overflow().SumW()
			default:
				weight_s_i1 = smearHs.Mlblb.Value(binx1)
			}

			switch binx2 {
			case hbook.UnderflowBin1D:
				weight_s_i2 = smearHs.Mlblb.Binning.Underflow().SumW()
			case hbook.OverflowBin1D:
				weight_s_i2 = smearHs.Mlblb.Binning.Overflow().SumW()
			default:
				weight_s_i2 = smearHs.Mlblb.Value(binx2)
			}

			Vec_Top[i_jets][0] = Vec_Top[i_jets][0] + weight_s_i1*weight_s_i2*Top_reco_px
			Vec_Top[i_jets][1] = Vec_Top[i_jets][1] + weight_s_i1*weight_s_i2*Top_reco_py
			Vec_Top[i_jets][2] = Vec_Top[i_jets][2] + weight_s_i1*weight_s_i2*Top_reco_pz

			Vec_Topbar[i_jets][0] = Vec_Topbar[i_jets][0] + weight_s_i1*weight_s_i2*Topbar_reco_px
			Vec_Topbar[i_jets][1] = Vec_Topbar[i_jets][1] + weight_s_i1*weight_s_i2*Topbar_reco_py
			Vec_Topbar[i_jets][2] = Vec_Topbar[i_jets][2] + weight_s_i1*weight_s_i2*Topbar_reco_pz

			weight_s_sum += weight_s_i1 * weight_s_i2

			i_analyzed_event = 1

			if debug {
				log.Printf("   weight            : %g", weight_s_sum)
				log.Printf("   (px, py, pz)[t]   : %g, %g, %g", Top_reco_px, Top_reco_py, Top_reco_pz)
				log.Printf("   (px, py, pz)[tbar]: %g, %g, %g", Topbar_reco_px, Topbar_reco_py, Topbar_reco_pz)
			}
		}

		if debug {
			log.Printf(" Number of iteration with solutions : %d", nIterations)
		}

		weights_com = append(weights_com, weight_s_sum)
		//		if i_jets == 0 {
		//			lep0 = lep
		//			lepbar0 = lepbar
		//		}
		//		if i_jets == 1 {
		//			lep1 = lep
		//			lepbar1 = lepbar
		//		}
	}

	// check whether we found a useful solution
	if i_analyzed_event != 1 || len(weights_com) == 0 {
		return nil
	}

	if debug {
		log.Printf(" Weight sum of jet combinatorics (0, 1): %, %g", weights_com[0], weights_com[1])
	}

	var (
		top_p_sum    r3.Vec
		topbar_p_sum r3.Vec

		i_weight_sel = 0
	)

	if weights_com[0] > weights_com[1] {
		if weights_com[0] > 0. && Vec_Top[0][0] != 0 {

			top_p_sum.X = Vec_Top[0][0] / weights_com[0]
			top_p_sum.Y = Vec_Top[0][1] / weights_com[0]
			top_p_sum.Z = Vec_Top[0][2] / weights_com[0]

			topbar_p_sum.X = Vec_Topbar[0][0] / weights_com[0]
			topbar_p_sum.Y = Vec_Topbar[0][1] / weights_com[0]
			topbar_p_sum.Z = Vec_Topbar[0][2] / weights_com[0]

			//	lplus = lepbar0
			//	lminus = lep0

			i_weight_sel = 1

		}

	} else {
		if weights_com[1] > 0. && Vec_Top[1][0] != 0 {

			top_p_sum.X = Vec_Top[1][0] / weights_com[1]
			top_p_sum.Y = Vec_Top[1][1] / weights_com[1]
			top_p_sum.Z = Vec_Top[1][2] / weights_com[1]

			topbar_p_sum.X = Vec_Topbar[1][0] / weights_com[1]
			topbar_p_sum.Y = Vec_Topbar[1][1] / weights_com[1]
			topbar_p_sum.Z = Vec_Topbar[1][2] / weights_com[1]

			//	lplus = lepbar1
			//	lminus = lep1

			i_weight_sel = 1

		}
	}

	if i_weight_sel == 0 {
		return nil
	}

	var (
		Top_fin_M    = m_top
		Topbar_fin_M = m_topbar

		Top_fin_E    = math.Sqrt(r3.Norm2(top_p_sum) + Top_fin_M*Top_fin_M)
		Topbar_fin_E = math.Sqrt(r3.Norm2(topbar_p_sum) + Topbar_fin_M*Topbar_fin_M)

		Top_fin    = fmom.NewPxPyPzE(top_p_sum.X, top_p_sum.Y, top_p_sum.Z, Top_fin_E)
		Topbar_fin = fmom.NewPxPyPzE(topbar_p_sum.X, topbar_p_sum.Y, topbar_p_sum.Z, Topbar_fin_E)
	)

	return []fmom.PxPyPzE{Top_fin, Topbar_fin}
}

func iabs(i int) int {
	if i < 0 {
		return -i
	}
	return i
}

func newPtEtaPhiM(pt, eta, phi, m float64) fmom.PxPyPzE {
	var (
		p = fmom.NewPtEtaPhiM(pt, eta, phi, m)
		o fmom.PxPyPzE
	)
	o.Set(&p)
	return o
}

func newPtEtaPhiE(pt, eta, phi, e float64) fmom.PxPyPzE {
	pt = math.Abs(pt)
	sin, cos := math.Sincos(phi)
	return fmom.NewPxPyPzE(
		pt*cos, pt*sin, pt*math.Sinh(eta),
		e,
	)
}
