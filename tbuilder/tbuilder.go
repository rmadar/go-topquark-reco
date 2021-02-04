package tbuilder

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

// Particle masses (in GeV)
const (
	mW      = 80.379
	mWbar   = 80.379
	mb      = 4.18
	mbbar   = 4.18
	mNu     = 0.0
	mNubar  = 0.0
	mTop    = 172.5
	mTopbar = 172.5
	mLep    = 0.0
	mLepbar = 0.0
)

type builderMode byte

const (
	undefined builderMode = iota
	sonnMode
	ellMode
	allMode
)

// TopBuilder reconstructs ttbar pairs in dilepton final state.
// Two methods can be used:
//  * the Sonnenschein method (https://arxiv.org/abs/hep-ph/0603011)
//  * the ellipses method (https://arxiv.org/abs/1305.1878)
// Each of these method is implemented along a smearing of
// object kinematics in order to maximize the reconstruction
// efficiency.
type TopBuilder struct {
	smearer       *smearingHistos
	rnd           *rand.Rand
	smearN        int
	smearLepPt    bool
	smearLepTheta bool
	smearLepAzimu bool
	smearJetPt    bool
	smearJetTheta bool
	smearJetAzimu bool
	mode          builderMode
	debug         bool
}

// New returns a new TopBuilder object from the path to a ROOT file holding
// histograms used to smear object kinematics, and set of options.
func New(fname string, opts ...Option) (*TopBuilder, error) {

	// Configuration of the reconstruction
	cfg := newConfig()
	for _, opt := range opts {
		opt(cfg)
	}

	// Get smearing histograms
	sh, err := newSmearingHistos(fname, cfg.rndseed)
	if err != nil {
		return nil, fmt.Errorf("could not create smearing histograms: %w", err)
	}

	// Activate all smearing if smearAll is enabled
	if cfg.smearAll {
		cfg.smearLepPt = true
		cfg.smearLepTheta = true
		cfg.smearLepAzimu = true
		cfg.smearJetPt = true
		cfg.smearJetTheta = true
		cfg.smearJetAzimu = true
	}

	// Create the TopBuilder object
	tb := &TopBuilder{
		smearer:       sh,
		rnd:           rand.New(rand.NewSource(cfg.rndseed)),
		smearLepPt:    cfg.smearLepPt,
		smearLepTheta: cfg.smearLepTheta,
		smearLepAzimu: cfg.smearLepAzimu,
		smearJetPt:    cfg.smearJetPt,
		smearJetTheta: cfg.smearJetTheta,
		smearJetAzimu: cfg.smearJetAzimu,
		smearN:        cfg.smearN,
		debug:         cfg.debug,
		mode:          undefined,
	}

	// Set the number of iteration to 1.0 if smearing is disabeld
	noLepSmear := !tb.smearLepPt && !tb.smearLepTheta && !tb.smearLepAzimu
	noJetSmear := !tb.smearJetPt && !tb.smearJetTheta && !tb.smearJetAzimu
	if noLepSmear && noJetSmear {
		tb.smearN = 1
	}

	if tb.debug {
		fmt.Printf("\n\n")
		log.Printf("Top Reconstruction Configuration:\n")
		log.Printf("  - Smearing file      : %v\n", fname)
		log.Printf("  - Random number seed : %v\n", cfg.rndseed)
		log.Printf("  - Smearing iterations: %v\n", cfg.smearN)
		log.Printf("  - Smearing lep pt    : %v\n", cfg.smearLepPt)
		log.Printf("  - Smearing lep theta : %v\n", cfg.smearLepTheta)
		log.Printf("  - Smearing lep azimu : %v\n", cfg.smearLepAzimu)
		log.Printf("  - Smearing jet pt    : %v\n", cfg.smearJetPt)
		log.Printf("  - Smearing jet theta : %v\n", cfg.smearJetTheta)
		log.Printf("  - Smearing jet azimu : %v\n\n\n", cfg.smearJetAzimu)
		/*
			switch tb.mode {
			case sonnMode:
				log.Printf("  - Mode :               Sonnenschein\n\n\n")
			case ellMode:
				log.Printf("  - Mode :               Ellipsis\n\n\n")
			case allMode:
				log.Printf("  - Mode :               All methods\n\n\n")
			default:
				log.Printf("invalid builder mode %d", tb.mode)
				return nil, fmt.Errorf("invalid builder mode (v=%d)", tb.mode)
			}
		*/
	}

	// Return the reconstruction object
	return tb, nil
}

// smearingHistos holds the histograms used to smear 4-vectors.
type smearingHistos struct {
	PtLepEE     hbook.Rand1D
	PtLepMu     hbook.Rand1D
	ThetaLepEE  hbook.Rand1D
	ThetaLepMu  hbook.Rand1D
	PtJet       hbook.Rand1D
	ThetaJet    hbook.Rand1D
	ThetaJetBar hbook.Rand1D
	Mlblb       *hbook.H1D
}

// NewSmearingHistos returns a set of smearing histograms from the provided
// ROOT file name and a seed for the histogram distributions.
func newSmearingHistos(fname string, seed uint64) (*smearingHistos, error) {
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

	sh := &smearingHistos{
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

// Reco performs the reconstruction and return 2 4-momentum
// (top and anti-top) as well as a status of the reconstruction.
// Status = N (0) means the reconstruction (didn't) work with N succesful iterations.
// Four-momemtum and missing Et components are in GeV.
// rdnNumbers are tuples of the 12 numbers needed per smearing iteration.
// If rndNumbers is not empty, smearing is automatically activated
// for this function call for all kinematic variables. When several
// tuples are passed, the  number of smearing iterations is set
// to be number of tuples, ie len(rdnNumbers). The ordering of
// random numbers is the following:
//   rdn[12] = [12]float64{
//     scaleLep , scaleLepBar ,
//     thetaLep , thetaLepBar ,
//     azimuLep , azimuLepBar ,
//     scaleJet0, scaleJet0Bar,
//     thetaJet0, thetaJet0Bar,
//     azimuJet0, azimuJet0Bar,
//     scaleJet1, scaleJet1Bar,
//     thetaJet1, thetaJet1Bar,
//     azimuJet1, azimuJet1Bar,
//   }
func (tb *TopBuilder) reconstruct(
	lepTLV, lepbarTLV fmom.PxPyPzE, pdgIDLep, pdgIDLepBar int,
	jetTLV, jetbarTLV fmom.PxPyPzE, isbJet, isbJetbar bool,
	emissx, emissy float64, rdnNumbers ...[12]float64) ([2]fmom.PxPyPzE, [2]fmom.PxPyPzE, [2]int) {

	// Configuration variables
	var (
		smearHs       = tb.smearer
		rnd           = tb.rnd
		usrRndN       = false
		smearLepPt    = tb.smearLepPt
		smearLepTheta = tb.smearLepTheta
		smearLepAzimu = tb.smearLepAzimu
		smearJetPt    = tb.smearJetPt
		smearJetTheta = tb.smearJetTheta
		smearJetAzimu = tb.smearJetAzimu
		smearN        = tb.smearN
		debug         = tb.debug
	)

	// Overwrite configuration for user-defined rnd numbers.
	// In that case, all smearing kinemarics is applied.
	if len(rdnNumbers) > 0 {
		usrRndN = true
		smearN = len(rdnNumbers)
	}

	// Ouptut of the algorithm
	var (
		tFinal    [2]fmom.PxPyPzE
		tbarFinal [2]fmom.PxPyPzE
		status    [2]int
	)

	// run the reconstruction only if the two jets are btags.
	if !isbJet || !isbJetbar {
		return tFinal, tbarFinal, status
	}

	var (
		lepIsE     = iabs(pdgIDLep) == 11
		lepIsMu    = iabs(pdgIDLep) == 13
		lepBarIsE  = iabs(pdgIDLepBar) == 11
		lepBarIsMu = iabs(pdgIDLepBar) == 13
	)

	var (
		Vec_Top     [2][2]r3.Vec
		Vec_Topbar  [2][2]r3.Vec
		weights_com [2][2]float64
		nIterations [2][2]int
	)

	// loop over jets
	for i_jets := 0; i_jets < 2; i_jets++ {

		if debug {
			log.Printf(" Jets combination %d:", i_jets)
		}

		// Number of iteration (for each reco method)
		nIterations[i_jets] = [2]int{0, 0}

		// Slices for the smeared kimenatics
		var (
			l   = make([]fmom.PxPyPzE, smearN)
			lb  = make([]fmom.PxPyPzE, smearN)
			j   = make([]fmom.PxPyPzE, smearN)
			jb  = make([]fmom.PxPyPzE, smearN)
			etx = make([]float64, smearN)
			ety = make([]float64, smearN)
		)

		// Starting smearing loop
		for i_smear := 0; i_smear < smearN; i_smear++ {

			// All four momenta for smearing
			var (
				jet, jetbar                   fmom.PxPyPzE
				lep, lepbar                   fmom.PxPyPzE
				lep_pt_smear, lepbar_pt_smear fmom.PxPyPzE
				lep_nosmear, lepbar_nosmear   fmom.PxPyPzE
				jet_pt_smear, jetbar_pt_smear fmom.PxPyPzE
				jet_nosmear, jetbar_nosmear   fmom.PxPyPzE
			)

			// All random variables to smear variables
			var (
				smear_scale_lep_0 = 1.0
				smear_scale_lep_1 = 1.0
				smear_theta_lep_0 = 0.0
				smear_theta_lep_1 = 0.0
				smear_azimu_lep_0 = 0.0
				smear_azimu_lep_1 = 0.0
				smear_scale_jet_0 = 1.0
				smear_scale_jet_1 = 1.0
				smear_theta_jet_0 = 0.0
				smear_theta_jet_1 = 0.0
				smear_azimu_jet_0 = 0.0
				smear_azimu_jet_1 = 0.0
			)

			if smearLepPt {
				switch {
				case lepIsE:
					smear_scale_lep_0 = smearHs.PtLepEE.Rand()
				case lepIsMu:
					smear_scale_lep_0 = smearHs.PtLepMu.Rand()
				default:
					smear_scale_lep_0 = 1
				}

				switch {
				case lepBarIsE:
					smear_scale_lep_1 = smearHs.PtLepEE.Rand()
				case lepBarIsMu:
					smear_scale_lep_1 = smearHs.PtLepMu.Rand()
				default:
					smear_scale_lep_1 = 1
				}
			}

			// User-defined random numbers
			if usrRndN {
				smear_scale_lep_0 = rdnNumbers[i_smear][0]
				smear_scale_lep_1 = rdnNumbers[i_smear][1]
			}

			lep_pt_smear.SetPtEtaPhiM(
				lepTLV.Pt()*smear_scale_lep_0,
				lepTLV.Eta(),
				lepTLV.Phi(),
				mLepbar,
			)

			lepbar_pt_smear.SetPtEtaPhiM(
				lepbarTLV.Pt()*smear_scale_lep_1,
				lepbarTLV.Eta(),
				lepbarTLV.Phi(),
				mLepbar,
			)

			lep_nosmear = lepTLV
			lepbar_nosmear = lepbarTLV

			// let's define the axis transverse to lepton momentum
			// for theta smearing accounting.
			var (
				lep_v_pt_1    = fmom.VecOf(&lep_pt_smear)
				lepbar_v_pt_1 = fmom.VecOf(&lepbar_pt_smear)
				beam_axis1    = r3.Vec{Z: 1}

				lep_axis   = lep_v_pt_1
				lep_axis_t = beam_axis1.Cross(lep_axis)

				lepbar_axis   = lepbar_v_pt_1
				lepbar_axis_t = beam_axis1.Cross(lepbar_axis)

				lep_v_pt    = lep_v_pt_1
				lepbar_v_pt = lepbar_v_pt_1
			)

			if smearLepTheta {
				switch {
				case lepIsE:
					smear_theta_lep_0 = smearHs.ThetaLepEE.Rand()
				case lepIsMu:
					smear_theta_lep_0 = smearHs.ThetaLepMu.Rand()
				}
				switch {
				case lepBarIsE:
					smear_theta_lep_1 = smearHs.ThetaLepEE.Rand()
				case lepBarIsMu:
					smear_theta_lep_1 = smearHs.ThetaLepMu.Rand()
				}
			}

			if smearLepAzimu {
				smear_azimu_lep_0 = rnd.Float64() * 2 * math.Pi
				smear_azimu_lep_1 = rnd.Float64() * 2 * math.Pi
			}

			// User-defined random numbers
			if usrRndN {
				smear_theta_lep_0 = rdnNumbers[i_smear][2]
				smear_theta_lep_1 = rdnNumbers[i_smear][3]
				smear_azimu_lep_0 = rdnNumbers[i_smear][4]
				smear_azimu_lep_1 = rdnNumbers[i_smear][5]
			}

			lep_v_pt = lep_v_pt.Rotate(smear_theta_lep_0, lep_axis_t)
			lepbar_v_pt = lepbar_v_pt.Rotate(smear_theta_lep_1, lepbar_axis_t)

			lep_v_pt = lep_v_pt.Rotate(smear_azimu_lep_0, lep_v_pt_1)
			lepbar_v_pt = lepbar_v_pt.Rotate(smear_azimu_lep_1, lepbar_v_pt_1)

			lep = fmom.NewPxPyPzE(
				lep_v_pt.X, lep_v_pt.Y, lep_v_pt.Z,
				lep_pt_smear.E(),
			)
			lepbar = fmom.NewPxPyPzE(
				lepbar_v_pt.X, lepbar_v_pt.Y, lepbar_v_pt.Z,
				lepbar_pt_smear.E(),
			)

			// Jet pT smearing
			if smearJetPt {
				smear_scale_jet_0 = smearHs.PtJet.Rand()
				smear_scale_jet_1 = smearHs.PtJet.Rand()
			}

			// Jet polar angle smearing
			if smearJetTheta {
				smear_theta_jet_0 = smearHs.ThetaJet.Rand()
				smear_theta_jet_1 = smearHs.ThetaJetBar.Rand()
			}

			// Jet Azimu angle smearing
			if smearJetAzimu {
				smear_azimu_jet_0 = rnd.Float64() * 2 * math.Pi
				smear_azimu_jet_1 = rnd.Float64() * 2 * math.Pi
			}

			// User-defined random numbers
			if usrRndN {
				smear_scale_jet_0 = rdnNumbers[i_smear][6]
				smear_scale_jet_1 = rdnNumbers[i_smear][7]
				smear_theta_jet_0 = rdnNumbers[i_smear][8]
				smear_theta_jet_1 = rdnNumbers[i_smear][9]
				smear_azimu_jet_0 = rdnNumbers[i_smear][10]
				smear_azimu_jet_1 = rdnNumbers[i_smear][11]
			}

			// To keep track of which random number for each jets
			// (mostly for xcheck when rnd numbers are passed by
			// the user).
			var (
				thetaJet    = 0.
				thetaJetBar = 0.
				azimuJet    = 0.
				azimuJetBar = 0.
			)

			switch i_jets {
			case 0:
				jet_pt_smear.SetPtEtaPhiE(
					jetTLV.Pt()*smear_scale_jet_0,
					jetTLV.Eta(),
					jetTLV.Phi(),
					jetTLV.E(),
				)
				jetbar_pt_smear.SetPtEtaPhiE(
					jetbarTLV.Pt()*smear_scale_jet_1,
					jetbarTLV.Eta(),
					jetbarTLV.Phi(),
					jetbarTLV.E(),
				)
				jet_nosmear = jetTLV
				jetbar_nosmear = jetbarTLV

				thetaJet = smear_theta_jet_0
				thetaJetBar = smear_theta_jet_1
				azimuJet = smear_azimu_jet_0
				azimuJetBar = smear_azimu_jet_1

			case 1:
				jetbar_pt_smear.SetPtEtaPhiE(
					jetTLV.Pt()*smear_scale_jet_0,
					jetTLV.Eta(),
					jetTLV.Phi(),
					jetTLV.E(),
				)

				jet_pt_smear.SetPtEtaPhiE(
					jetbarTLV.Pt()*smear_scale_jet_1,
					jetbarTLV.Eta(),
					jetbarTLV.Phi(),
					jetbarTLV.E(),
				)

				jetbar_nosmear = jetTLV
				jet_nosmear = jetbarTLV

				thetaJetBar = smear_theta_jet_0
				thetaJet = smear_theta_jet_1
				azimuJetBar = smear_azimu_jet_0
				azimuJet = smear_azimu_jet_1
			}

			// let's define the transverse axis to the jet momentum
			// for theta smearing accounting.
			var (
				jet_v_pt_1    = fmom.VecOf(&jet_pt_smear)
				jetbar_v_pt_1 = fmom.VecOf(&jetbar_pt_smear)

				jet_axis   = jet_v_pt_1
				jet_axis_t = beam_axis1.Cross(jet_axis)

				jetbar_axis   = jetbar_v_pt_1
				jetbar_axis_t = beam_axis1.Cross(jetbar_axis)

				jet_v_pt    = jet_v_pt_1
				jetbar_v_pt = jetbar_v_pt_1
			)

			jet_v_pt = jet_v_pt.Rotate(thetaJet, jet_axis_t)
			jetbar_v_pt = jetbar_v_pt.Rotate(thetaJetBar, jetbar_axis_t)

			jet_v_pt = jet_v_pt.Rotate(azimuJet, jet_v_pt_1)
			jetbar_v_pt = jetbar_v_pt.Rotate(azimuJetBar, jetbar_v_pt_1)

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
				log.Printf("   Lep    smearing (scale, theta, azimu): %7.5f, %7.5f %7.5f",
					smear_scale_lep_0, smear_theta_lep_0, smear_azimu_lep_0)
				log.Printf("   Lepbar smearing (scale, theta, azimu): %7.5f, %7.5f %7.5f",
					smear_scale_lep_1, smear_theta_lep_1, smear_azimu_lep_1)
				log.Printf("   Jet0   smearing (scale, theta, azimu): %7.5f, %7.5f %7.5f",
					smear_scale_jet_0, smear_theta_jet_0, smear_azimu_jet_0)
				log.Printf("   Jet1   smearing (scale, theta, azimu): %7.5f, %7.5f %7.5f",
					smear_scale_jet_1, smear_theta_jet_1, smear_azimu_jet_1)
				log.Printf("   Lep    Pt, Px, Py, Pz : %5.2f, %5.2f, %5.2f, %5.2f",
					lep_nosmear.Pt(),
					lep_nosmear.Px(),
					lep_nosmear.Py(),
					lep_nosmear.Pz(),
				)
				log.Printf("                   after : %5.3f, %5.3f, %5.3f, %5.3f",
					lep.Pt(),
					lep.Px(),
					lep.Py(),
					lep.Pz(),
				)
				log.Printf("   LepBar Pt, Px, Py, Pz : %5.3f, %5.3f, %5.3f, %5.3f",
					lepbar_nosmear.Pt(),
					lepbar_nosmear.Px(),
					lepbar_nosmear.Py(),
					lepbar_nosmear.Pz(),
				)
				log.Printf("                   after : %5.3f, %5.3f, %5.3f, %5.3f",
					lepbar.Pt(),
					lepbar.Px(),
					lepbar.Py(),
					lepbar.Pz(),
				)
				log.Printf("   Jet0   Pt, Px, Py, Pz : %5.3f, %5.3f, %5.3f, %5.3f",
					jet_nosmear.Pt(),
					jet_nosmear.Px(),
					jet_nosmear.Py(),
					jet_nosmear.Pz(),
				)
				log.Printf("                   after : %5.3f, %5.3f, %5.3f, %5.3f",
					jet.Pt(),
					jet.Px(),
					jet.Py(),
					jet.Pz(),
				)
				log.Printf("   Jet1   Pt, Px, Py, Pz : %5.3f, %5.3f, %5.3f, %5.3f",
					jetbar_nosmear.Pt(),
					jetbar_nosmear.Px(),
					jetbar_nosmear.Py(),
					jetbar_nosmear.Pz(),
				)
				log.Printf("                   after : %5.3f, %5.3f, %5.3f, %5.3f",
					jetbar.Pt(),
					jetbar.Px(),
					jetbar.Py(),
					jetbar.Pz(),
				)
				log.Printf("   MET        Pt, Px, Py : %5.3f, %5.3f, %5.3f",
					math.Hypot(Emiss_x_nosmear, Emiss_y_nosmear),
					Emiss_x_nosmear,
					Emiss_y_nosmear,
				)
				log.Printf("                   after : %5.3f, %5.3f, %5.3f",
					math.Hypot(Emiss_x_smear, Emiss_y_smear),
					Emiss_x_smear,
					Emiss_y_smear,
				)
			}

			// Fill the slices with smeared kinematics quantites
			l[i_smear] = lep
			lb[i_smear] = lepbar
			j[i_smear] = jet
			jb[i_smear] = jetbar
			etx[i_smear] = Emiss_x_smear
			ety[i_smear] = Emiss_y_smear
		}

		// Sonenschein: run the reconstruction over all smearing iterations and combined them.
		if tb.mode == sonnMode || tb.mode == allMode {
			p3t, p3tbar, w, nIter := recoCombineIters(l, lb, j, jb, etx, ety, sonnMode, tb.smearer, tb.debug)
			recoMode := 0
			Vec_Top[i_jets][recoMode] = p3t
			Vec_Topbar[i_jets][recoMode] = p3tbar
			nIterations[i_jets][recoMode] = nIter
			weights_com[recoMode][i_jets] = w
		}

		// Ellipse: run the reconstruction over all smearing iterations and combined them.
		if tb.mode == ellMode || tb.mode == allMode {
			p3t, p3tbar, w, nIter := recoCombineIters(l, lb, j, jb, etx, ety, ellMode, tb.smearer, tb.debug)
			recoMode := 1
			Vec_Top[i_jets][recoMode] = p3t
			Vec_Topbar[i_jets][recoMode] = p3tbar
			nIterations[i_jets][recoMode] = nIter
			weights_com[recoMode][i_jets] = w
		}

		if debug {
			log.Printf("   number of iteration with solutions : %d", nIterations[i_jets])
		}
	}

	// Check whether we found a useful solution,
	// depending of which reco were ran.
	zeros := [2]float64{0.0, 0.0}
	sonFails := weights_com[0] == zeros && tb.mode == sonnMode
	ellFails := weights_com[1] == zeros && tb.mode == ellMode
	allFails := weights_com[0] == zeros && weights_com[1] == zeros && tb.mode == allMode
	if sonFails || ellFails || allFails {
		return tFinal, tbarFinal, status
	}

	if debug {
		log.Printf(" Weight sum of jet combinatorics (0, 1): %5.2e, %5.2e", weights_com[0], weights_com[1])
	}

	// Jet combination selection, for each reconstruction.
	methodName := [2]string{"Sonnenschein", "Ellipsis"}
	for iReco := 0; iReco < 2; iReco++ {

		// Just skip, if there is no weights for the current method
		if weights_com[iReco] == zeros {
			continue
		}

		// Get the best jet combinatorics index
		jetComb := -1
		wJet1, wJet2 := weights_com[iReco][0], weights_com[iReco][1]
		if wJet1 > wJet2 {
			jetComb = 0
		} else if wJet2 > 0. && Vec_Top[1][iReco].X != 0 {
			jetComb = 1
		} else {
			continue
		}

		// Get the final top 4-momenta for the current method
		pt, ptb := Vec_Top[jetComb][iReco], Vec_Topbar[jetComb][iReco]
		Top_fin_E := math.Sqrt(r3.Norm2(pt) + mTop*mTop)
		Topbar_fin_E := math.Sqrt(r3.Norm2(ptb) + mTopbar*mTopbar)
		tFinal[iReco] = fmom.NewPxPyPzE(pt.X, pt.Y, pt.Z, Top_fin_E)
		tbarFinal[iReco] = fmom.NewPxPyPzE(ptb.X, ptb.Y, ptb.Z, Topbar_fin_E)

		// Status of the reconstruction
		if !isBad(tFinal[iReco]) && !isBad(tbarFinal[iReco]) {
			status[iReco] = nIterations[jetComb][iReco]
		}

		if debug {
			log.Printf(" %v method:\n", methodName[iReco])
			log.Printf("  Final P[t]   : %v, M=%v", tFinal[iReco], tFinal[iReco].M())
			log.Printf("  Final P[tbar]: %v, M=%v", tbarFinal[iReco], tbarFinal[iReco].M())
		}
	}

	// Return the results
	return tFinal, tbarFinal, status
}

// Helper function reconstructing & combining various
// smearing iterations, by doing the weighted mean.
// The weight is obtained from M(l,b) distribution.
func recoCombineIters(
	l, lb, j, jb []fmom.PxPyPzE, etx, ety []float64,
	mode builderMode, smearHs *smearingHistos, debug bool) (r3.Vec, r3.Vec, float64, int) {

	// Define output container for the return
	var (
		weight      = 0.0
		nIter       = 0
		p3t, p3tbar r3.Vec
	)

	// Define the reco method
	var reco func(
		fmom.PxPyPzE, fmom.PxPyPzE,
		fmom.PxPyPzE, fmom.PxPyPzE,
		float64, float64, bool) (r3.Vec, r3.Vec, bool)
	switch mode {
	case sonnMode:
		reco = sonnenschein
	case ellMode:
		reco = ellipsis
	default:
		panic("FIX THIS ERROR MESSAGE")
	}

	// Loop over iterations
	for i := range l {

		// The reco
		t, tbar, OK := reco(l[i], lb[i], j[i], jb[i], etx[i], ety[i], debug)

		// Get the weight
		var (
			mlbarb fmom.PxPyPzE
			mlbbar fmom.PxPyPzE
		)
		mlbarb.Set(fmom.Add(&lb[i], &j[i]))
		mlbbar.Set(fmom.Add(&l[i], &jb[i]))

		var (
			binx1 = hbook.Bin1Ds(smearHs.Mlblb.Binning.Bins).IndexOf(mlbarb.M())
			binx2 = hbook.Bin1Ds(smearHs.Mlblb.Binning.Bins).IndexOf(mlbbar.M())
			w1    = 0.0
			w2    = 0.0
		)

		// Only if the reconstruction is successful:
		//  - increment the number of good iteration
		//  - compute the weight associated to this jet combination
		if OK {
			nIter += 1

			switch binx1 {
			case hbook.UnderflowBin1D:
				w1 = smearHs.Mlblb.Binning.Underflow().SumW()
			case hbook.OverflowBin1D:
				w1 = smearHs.Mlblb.Binning.Overflow().SumW()
			default:
				w1 = smearHs.Mlblb.Value(binx1)
			}

			switch binx2 {
			case hbook.UnderflowBin1D:
				w2 = smearHs.Mlblb.Binning.Underflow().SumW()
			case hbook.OverflowBin1D:
				w2 = smearHs.Mlblb.Binning.Overflow().SumW()
			default:
				w2 = smearHs.Mlblb.Value(binx2)
			}
		}

		p3t = p3t.Add(t.Scale(w1 * w2))
		p3tbar = p3tbar.Add(tbar.Scale(w1 * w2))
		weight += w1 * w2

		if debug {
			log.Printf("   (px, py, pz)[t]   : %3.2f, %3.2f, %3.2f", p3t.X, p3t.Y, p3t.Z)
			log.Printf("   (px, py, pz)[tbar]: %3.2f, %3.2f, %3.2f", p3tbar.X, p3tbar.Y, p3tbar.Z)
			log.Printf("   weight[t, tbar]   : %5.3e, %5.3e", w1, w2)
			log.Printf("   weight sum        : %5.3e", weight)
		}

	}

	// Divide by the total weight to get the proper average
	p3t = p3t.Scale(1. / weight)
	p3tbar = p3tbar.Scale(1. / weight)

	// Return the result
	return p3t, p3tbar, weight, nIter
}

// SonnReco performs the Sonnenschein reconstruction and return 2 4-momentum
// (top and anti-top) as well as a status of the reconstruction.
// Status = N means the reconstruction worked for N smearing iterations.
// Therefore, Status = 0 means the reconstruction has failed.
// Four-momemtum and missing Et components must be given in GeV.
// rdnNumbers are tuples of the 12 numbers needed per smearing iteration.
// If rndNumbers is not empty, smearing is automatically activated
// for this function call for all kinematic variables. When several
// tuples are passed, the  number of smearing iterations is set
// to be number of tuples, ie len(rdnNumbers). The ordering of
// random numbers is the following:
//   rdn[12] = [12]float64{
//     scaleLep , scaleLepBar ,
//     thetaLep , thetaLepBar ,
//     azimuLep , azimuLepBar ,
//     scaleJet0, scaleJet0Bar,
//     thetaJet0, thetaJet0Bar,
//     azimuJet0, azimuJet0Bar,
//     scaleJet1, scaleJet1Bar,
//     thetaJet1, thetaJet1Bar,
//     azimuJet1, azimuJet1Bar,
//   }
func (tb *TopBuilder) SonnReco(
	lepTLV, lepbarTLV fmom.PxPyPzE, pdgIDLep, pdgIDLepBar int,
	jetTLV, jetbarTLV fmom.PxPyPzE, isbJet, isbJetbar bool,
	emissx, emissy float64, rdnNumbers ...[12]float64) (fmom.PxPyPzE, fmom.PxPyPzE, int) {

	// Set mode to Sonnenschein reco
	tb.mode = sonnMode

	// Run the reco
	top, antiop, status := tb.reconstruct(
		lepTLV, lepbarTLV, pdgIDLep,
		pdgIDLepBar, jetTLV, jetbarTLV,
		isbJet, isbJetbar, emissx, emissy,
		rdnNumbers...,
	)

	// Return the Sonnenschein result only
	return top[0], antiop[0], status[0]
}

// ElliReco performs the Sonnenschein reconstruction and return 2 4-momentum
// (top and anti-top) as well as a status of the reconstruction.
// Status = N means the reconstruction worked for N smearing iterations.
// Therefore, Status = 0 means the reconstruction has failed.
// Four-momemtum and missing Et components must be given in GeV.
// rdnNumbers are tuples of the 12 numbers needed per smearing iteration.
// If rndNumbers is not empty, smearing is automatically activated
// for this function call for all kinematic variables. When several
// tuples are passed, the  number of smearing iterations is set
// to be number of tuples, ie len(rdnNumbers). The ordering of
// random numbers is the following:
//   rdn[12] = [12]float64{
//     scaleLep , scaleLepBar ,
//     thetaLep , thetaLepBar ,
//     azimuLep , azimuLepBar ,
//     scaleJet0, scaleJet0Bar,
//     thetaJet0, thetaJet0Bar,
//     azimuJet0, azimuJet0Bar,
//     scaleJet1, scaleJet1Bar,
//     thetaJet1, thetaJet1Bar,
//     azimuJet1, azimuJet1Bar,
//   }
func (tb *TopBuilder) ElliReco(
	lepTLV, lepbarTLV fmom.PxPyPzE, pdgIDLep, pdgIDLepBar int,
	jetTLV, jetbarTLV fmom.PxPyPzE, isbJet, isbJetbar bool,
	emissx, emissy float64, rdnNumbers ...[12]float64) (fmom.PxPyPzE, fmom.PxPyPzE, int) {

	// Set mode to Ellipse method
	tb.mode = ellMode

	// Run the reco
	top, antiop, status := tb.reconstruct(
		lepTLV, lepbarTLV, pdgIDLep,
		pdgIDLepBar, jetTLV, jetbarTLV,
		isbJet, isbJetbar, emissx, emissy,
		rdnNumbers...,
	)

	// Return the ellipse result only
	return top[1], antiop[1], status[1]
}

// AllReco performs the Sonnenschein and the Ellipse reconstruction
// and return 2 array of 4-momenta (top[2] and anti-top[2]) as well as
// an array of status[2] of the reconstruction. The array index value
// is 0 for the Sonnenschein method and 1 for the Ellipse method.
// Status[i] = N means the reconstruction i worked for N smearing iterations.
// Therefore, Status[i] = 0 means the reconstruction has failed.
// Four-momemtum and missing Et components must be given in GeV.
// rdnNumbers are tuples of the 12 numbers needed per smearing iteration.
// If rndNumbers is not empty, smearing is automatically activated
// for this function call for all kinematic variables. When several
// tuples are passed, the  number of smearing iterations is set
// to be number of tuples, ie len(rdnNumbers). The ordering of
// random numbers is the following:
//   rdn[12] = [12]float64{
//     scaleLep , scaleLepBar ,
//     thetaLep , thetaLepBar ,
//     azimuLep , azimuLepBar ,
//     scaleJet0, scaleJet0Bar,
//     thetaJet0, thetaJet0Bar,
//     azimuJet0, azimuJet0Bar,
//     scaleJet1, scaleJet1Bar,
//     thetaJet1, thetaJet1Bar,
//     azimuJet1, azimuJet1Bar,
//   }
func (tb *TopBuilder) AllReco(
	lepTLV, lepbarTLV fmom.PxPyPzE, pdgIDLep, pdgIDLepBar int,
	jetTLV, jetbarTLV fmom.PxPyPzE, isbJet, isbJetbar bool,
	emissx, emissy float64, rdnNumbers ...[12]float64) ([2]fmom.PxPyPzE, [2]fmom.PxPyPzE, [2]int) {

	// Set mode to all reco
	tb.mode = allMode

	// Return the both results
	return tb.reconstruct(
		lepTLV, lepbarTLV, pdgIDLep,
		pdgIDLepBar, jetTLV, jetbarTLV,
		isbJet, isbJetbar, emissx, emissy,
		rdnNumbers...,
	)
}

// Helper function returning the absolute value of an int.
func iabs(i int) int {
	if i < 0 {
		return -i
	}
	return i
}

// Helper function checking of the P4[top] makes sense.
func isBad(t fmom.PxPyPzE) bool {
	isDefault := t.Px() == 10000. && t.Py() == 10000. && t.Pz() == 10000.
	isEmpty := t.Px() == 0. && t.Py() == 0. && t.Pz() == 0.
	return isDefault || isEmpty
}

// Internal function performing the Sonnenschein reconsruction.
func sonnenschein(
	lep, lepbar, jet, jetbar fmom.PxPyPzE, etx, ety float64,
	debug bool) (r3.Vec, r3.Vec, bool) {

	var (
		jet_v    = fmom.VecOf(&jet)
		lep_v    = fmom.VecOf(&lep)
		jetbar_v = fmom.VecOf(&jetbar)
		lepbar_v = fmom.VecOf(&lepbar)

		p_nu_y       = -999.0
		p_nu_z       = -999.0
		p_nubar_x    = -999.0
		p_nubar_y    = -999.0
		p_nubar_z    = -999.0
		p_nu_x_close = 100000.0

		Top_calc, Topbar_calc fmom.PxPyPzE
		nu_calc, nubar_calc   fmom.PxPyPzE
		Top1, Top2            fmom.PxPyPzE
		nu, top, topbar       fmom.PxPyPzE
		mtt_val               []float64
	)

	a1 := (jet.E()+lepbar.E())*(mW*mW-mLepbar*mLepbar-mNu*mNu) -
		lepbar.E()*(mTop*mTop-mb*mb-mLepbar*mLepbar-mNu*mNu) +
		2*jet.E()*lepbar.E()*lepbar.E() - 2*lepbar.E()*(jet_v.Dot(lepbar_v))
	a2 := 2 * (jet.E()*lepbar.Px() - lepbar.E()*jet.Px())
	a3 := 2 * (jet.E()*lepbar.Py() - lepbar.E()*jet.Py())
	a4 := 2 * (jet.E()*lepbar.Pz() - lepbar.E()*jet.Pz())

	b1 := (jetbar.E()+lep.E())*(mWbar*mWbar-mLep*mLep-mNubar*mNubar) -
		lep.E()*(mTopbar*mTopbar-mbbar*mbbar-mLep*mLep-mNubar*mNubar) +
		2.*jetbar.E()*lep.E()*lep.E() - 2.*lep.E()*(jetbar_v.Dot(lep_v))
	b2 := 2. * (jetbar.E()*lep.Px() - lep.E()*jetbar.Px())
	b3 := 2. * (jetbar.E()*lep.Py() - lep.E()*jetbar.Py())
	b4 := 2. * (jetbar.E()*lep.Pz() - lep.E()*jetbar.Pz())

	c22 := (mW*mW-mLepbar*mLepbar-mNu*mNu)*(mW*mW-mLepbar*mLepbar-mNu*mNu) -
		4.*(lepbar.E()*lepbar.E()-lepbar.Pz()*lepbar.Pz())*(a1/a4)*(a1/a4) -
		4.*(mW*mW-mLepbar*mLepbar-mNu*mNu)*lepbar.Pz()*a1/a4
	c21 := 4.*(mW*mW-mLepbar*mLepbar-mNu*mNu)*(lepbar.Px()-lepbar.Pz()*a2/a4) -
		8.*(lepbar.E()*lepbar.E()-lepbar.Pz()*lepbar.Pz())*a1*a2/(a4*a4) -
		8*lepbar.Px()*lepbar.Pz()*a1/a4
	c20 := -4.*(lepbar.E()*lepbar.E()-lepbar.Px()*lepbar.Px()) -
		4.*(lepbar.E()*lepbar.E()-lepbar.Pz()*lepbar.Pz())*(a2/a4)*(a2/a4) -
		8.*lepbar.Px()*lepbar.Pz()*a2/a4

	c11 := 4.*(mW*mW-mLepbar*mLepbar-mNu*mNu)*(lepbar.Py()-lepbar.Pz()*a3/a4) -
		8.*(lepbar.E()*lepbar.E()-lepbar.Pz()*lepbar.Pz())*a1*a3/(a4*a4) -
		8.*lepbar.Py()*lepbar.Pz()*a1/a4
	c10 := -8.*(lepbar.E()*lepbar.E()-lepbar.Pz()*lepbar.Pz())*a2*a3/(a4*a4) +
		8.*lepbar.Px()*lepbar.Py() - 8.*lepbar.Px()*lepbar.Pz()*a3/a4 -
		8.*lepbar.Py()*lepbar.Pz()*a2/a4
	c00 := -4.*(lepbar.E()*lepbar.E()-lepbar.Py()*lepbar.Py()) -
		4.*(lepbar.E()*lepbar.E()-lepbar.Pz()*lepbar.Pz())*(a3/a4)*(a3/a4) -
		8.*lepbar.Py()*lepbar.Pz()*a3/a4

	if debug {
		log.Printf("   a1, a2, a3, a4 = %5.3f, %5.3f, %5.3f, %5.3f", a1, a2, a3, a4)
		log.Printf("   b1, b2, b3, b4 = %5.3f, %5.3f, %5.3f, %5.3f", b1, b2, b3, b4)
		log.Printf("   c00, c10, c11  = %5.3f, %5.3f, %5.3f", c00, c10, c11)
		log.Printf("   c22, c21, c20  = %5.3f, %5.3f, %5.3f", c22, c21, c20)
	}

	// new norm
	n44 := a4 * a4
	c22 *= n44
	c21 *= n44
	c20 *= n44
	c11 *= n44
	c10 *= n44
	c00 *= n44

	dp22 := (mWbar*mWbar-mLep*mLep-mNubar*mNubar)*(mWbar*mWbar-mLep*mLep-mNubar*mNubar) -
		4.*(lep.E()*lep.E()-lep.Pz()*lep.Pz())*(b1/b4)*(b1/b4) -
		4.*(mWbar*mWbar-mLep*mLep-mNubar*mNubar)*lep.Pz()*b1/b4
	dp21 := 4.*(mWbar*mWbar-mLep*mLep-mNubar*mNubar)*(lep.Px()-lep.Pz()*b2/b4) -
		8.*(lep.E()*lep.E()-lep.Pz()*lep.Pz())*b1*b2/(b4*b4) - 8.*lep.Px()*lep.Pz()*b1/b4
	dp20 := -4.*(lep.E()*lep.E()-lep.Px()*lep.Px()) -
		4.*(lep.E()*lep.E()-lep.Pz()*lep.Pz())*(b2/b4)*(b2/b4) -
		8.*lep.Px()*lep.Pz()*b2/b4

	dp11 := 4*(mWbar*mWbar-mLep*mLep-mNubar*mNubar)*(lep.Py()-lep.Pz()*b3/b4) -
		8.*(lep.E()*lep.E()-lep.Pz()*lep.Pz())*b1*b3/(b4*b4) -
		8.*lep.Py()*lep.Pz()*b1/b4
	dp10 := -8.*(lep.E()*lep.E()-lep.Pz()*lep.Pz())*b2*b3/(b4*b4) +
		8.*lep.Px()*lep.Py() - 8.*lep.Px()*lep.Pz()*b3/b4 - 8.*lep.Py()*lep.Pz()*b2/b4
	dp00 := -4.*(lep.E()*lep.E()-lep.Py()*lep.Py()) -
		4.*(lep.E()*lep.E()-lep.Pz()*lep.Pz())*(b3/b4)*(b3/b4) -
		8.*lep.Py()*lep.Pz()*b3/b4

	d22 := dp22 +
		etx*etx*dp20 +
		ety*ety*dp00 +
		etx*ety*dp10 +
		etx*dp21 +
		ety*dp11
	d21 := -1.*dp21 - 2.*etx*dp20 - ety*dp10
	d20 := dp20
	d11 := -1.*dp11 - 2.*ety*dp00 - etx*dp10
	d10 := dp10
	d00 := dp00

	d22 = d22 * b4 * b4
	d21 = d21 * b4 * b4
	d20 = d20 * b4 * b4
	d11 = d11 * b4 * b4
	d10 = d10 * b4 * b4
	d00 = d00 * b4 * b4

	if debug {
		log.Printf("   d22, d21, d20  = %5.3f, %5.3f, %5.3f", d22, d21, d20)
		log.Printf("   d11, d10, d00  = %5.3f, %5.3f, %5.3f", d11, d10, d00)
	}

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

	if debug {
		log.Printf("   polynom coeff (h0, h1, h2, h3, h4) = %.3e, %.3e, %.3e, %.3e, %.3e\n", h0, h1, h2, h3, h4)
	}

	roots := make([]float64, 0, 4)
	const ε = 1e-15
	for _, z := range zs {
		if math.Abs(imag(z)) < ε {
			roots = append(roots, real(z))
		}
	}
	nSolutions := len(roots)
	if debug {
		log.Printf("   number of polynom roots: %d", nSolutions)
	}

	// Default values
	var (
		mtt_min = 100000.0

		Top_reco_px = 10000.0
		Top_reco_py = 10000.0
		Top_reco_pz = 10000.0

		Topbar_reco_px = 10000.0
		Topbar_reco_py = 10000.0
		Topbar_reco_pz = 10000.0
	)

	// If no solution, return default values.
	if nSolutions == 0 {
		var (
			t1 = r3.Vec{X: Top_reco_px, Y: Top_reco_py, Z: Top_reco_pz}
			t2 = r3.Vec{X: Topbar_reco_px, Y: Topbar_reco_py, Z: Topbar_reco_pz}
			ok bool
		)
		return t1, t2, ok
	}

	// Select solutions
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
		p_nubar_x = etx - p_nu_x_close
		p_nubar_y = ety - p_nu_y
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
			Top_reco_px = Top_calc.Px()
			Top_reco_py = Top_calc.Py()
			Top_reco_pz = Top_calc.Pz()
			Topbar_reco_px = Topbar_calc.Px()
			Topbar_reco_py = Topbar_calc.Py()
			Topbar_reco_pz = Topbar_calc.Pz()
		} // ---> end of solution selection
	}

	// Final results
	topP := r3.Vec{X: Top_reco_px, Y: Top_reco_py, Z: Top_reco_pz}
	topbarP := r3.Vec{X: Topbar_reco_px, Y: Topbar_reco_py, Z: Topbar_reco_pz}

	return topP, topbarP, nSolutions > 0
}
