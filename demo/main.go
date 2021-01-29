package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"strings"

	"go-hep.org/x/hep/fmom"
	"go-hep.org/x/hep/groot"
	"go-hep.org/x/hep/groot/rtree"

	"github.com/rmadar/go-topquark-reco/tbuilder"
)

func main() {

	// Command line arguments
	var (
		nEvts         = flag.Int("n", -1, "Number of processed events")
		debug         = flag.Bool("debug", false, "Print various numbers along the process")
		smearAll      = flag.Bool("smearAll", false, "Enable smearing of all quantities")
		smearN        = flag.Int("smearN", 100, "Number of smearing iterations")
		smearLepPt    = flag.Bool("smearLepPt", false, "Enable lepton pT smearing")
		smearLepTheta = flag.Bool("smearLepTheta", false, "Enable lepton polar angle smearing")
		smearLepAzimu = flag.Bool("smearLepAzimy", false, "Enable lepton azimuth angle smearing")
		smearJetPt    = flag.Bool("smearJetPt", false, "Enable jet pT smearing")
		smearJetTheta = flag.Bool("smearJetTheta", false, "Enable jet polar angle smearing")
		smearJetAzimu = flag.Bool("smearJetAzimy", false, "Enable jet azimuth angle smearing")
		reco          = flag.String("reco", "sonn", "Select reconstruction method (sonn[enschein]|ell[ipsis])")
	)
	flag.Parse()

	// Open the test ROOT file
	f, err := groot.Open("../testdata/data.root")
	if err != nil {
		log.Fatalf("could not open ROOT file: %+v", err)
	}
	defer f.Close()

	// Get the TTree
	o, err := f.Get("nominal")
	if err != nil {
		log.Fatalf("could not retrieve ROOT tree: %+v", err)
	}
	t := o.(rtree.Tree)

	// Get variables to read
	var (
		nBad      = 0
		evtNum    int64
		lepPt     []float32
		lepEta    []float32
		lepPhi    []float32
		lepPid    []int32
		jetPt     []float32
		jetE      []float32
		jetEta    []float32
		jetPhi    []float32
		jetMV2c10 []float32
		nBjets    int32
		metMet    float32
		metPhi    float32

		rvars = []rtree.ReadVar{
			{Name: "eventNumber", Value: &evtNum},
			{Name: "d_lep_pt", Value: &lepPt},
			{Name: "d_lep_eta", Value: &lepEta},
			{Name: "d_lep_phi", Value: &lepPhi},
			{Name: "d_lep_pid", Value: &lepPid},
			{Name: "d_jet_pt", Value: &jetPt},
			{Name: "d_jet_e", Value: &jetE},
			{Name: "d_jet_eta", Value: &jetEta},
			{Name: "d_jet_phi", Value: &jetPhi},
			{Name: "d_jet_mv2c10", Value: &jetMV2c10},
			{Name: "d_nbjet", Value: &nBjets},
			{Name: "d_met_met", Value: &metMet},
			{Name: "d_met_phi", Value: &metPhi},
		}
	)

	// Get the TTree reader
	r, err := rtree.NewReader(t, rvars, rtree.WithRange(0, int64(*nEvts)))
	if err != nil {
		log.Fatalf("could not create tree reader: %+v", err)
	}
	defer r.Close()

	var recoMeth func() tbuilder.Option
	switch mode := strings.ToLower(*reco); {
	case strings.HasPrefix(mode, "sonn"):
		recoMeth = tbuilder.WithSonnenschein
	case strings.HasPrefix(mode, "ell"):
		recoMeth = tbuilder.WithEllipsis
	default:
		log.Panicf("unknown reconstruction method %q", *reco)
	}

	// create a Sonnenschein reco algorithm.
	topBuilder, err := tbuilder.New("../testdata/smearingHistos.root",
		tbuilder.WithDebug(*debug),
		tbuilder.WithSmearN(*smearN),
		tbuilder.WithSmearAll(*smearAll),
		tbuilder.WithSmearLepPt(*smearLepPt),
		tbuilder.WithSmearLepTheta(*smearLepTheta),
		tbuilder.WithSmearLepAzimu(*smearLepAzimu),
		tbuilder.WithSmearJetPt(*smearJetPt),
		tbuilder.WithSmearJetTheta(*smearJetTheta),
		tbuilder.WithSmearJetAzimu(*smearJetAzimu),
		recoMeth(),
	)

	if err != nil {
		log.Fatalf("could not create the Sonnenschein reco algorithm: %+v", err)
	}

	// Event loop
	err = r.Read(func(ctx rtree.RCtx) error {

		// Prepare leptons four vectors
		var (
			lep0 fmom.PxPyPzE
			lep1 fmom.PxPyPzE
			jet0 fmom.PxPyPzE
			jet1 fmom.PxPyPzE

			lid0 int
			lid1 int
		)

		for i, id := range lepPid {
			if id > 0 {
				lid1 = int(id)
				lep1.SetPtEtaPhiM(float64(lepPt[i]), float64(lepEta[i]), float64(lepPhi[i]), 0.0)
			} else {
				lid0 = int(id)
				lep0.SetPtEtaPhiM(float64(lepPt[i]), float64(lepEta[i]), float64(lepPhi[i]), 0.0)
			}
		}

		// Prepare jet four vectors and b-tagg info based on the 2 leading jets
		jet0.SetPtEtaPhiE(float64(jetPt[0]), float64(jetEta[0]), float64(jetPhi[1]), float64(jetE[0]))
		jet1.SetPtEtaPhiE(float64(jetPt[1]), float64(jetEta[1]), float64(jetPhi[1]), float64(jetE[1]))
		var (
			j1b = jetMV2c10[0] > 0.691
			j2b = jetMV2c10[1] > 0.691
		)

		// Prepare missing transverse energy component
		sin, cos := math.Sincos(float64(metPhi))
		Etx := float64(metMet) * cos
		Ety := float64(metMet) * sin

		// Call the reconstruction of top and antitop
		t, tbar, nIterations := topBuilder.Reconstruct(
			lep0, lep1, lid0, lid1,
			jet0, jet1, j1b, j2b,
			Etx, Ety,
		)

		// Keep track of not reconstructed events
		if nIterations == 0 {
			nBad++
		}

		// Print some information
		fmt.Printf("Entry %d:\n", ctx.Entry)
		fmt.Printf("   - Evt number   %v\n", evtNum)
		fmt.Printf("   - N[b-jets]    %v\n", nBjets)
		fmt.Printf("   - b-tagging:   1st=%v, 2nd=%v\n", j1b, j2b)
		fmt.Printf("   - final state  %v\n", lepPid)
		fmt.Printf("   - P4[l]        %v\n", lep0)
		fmt.Printf("   - P4[lbar]     %v\n", lep1)
		fmt.Printf("   - P4[top]      %v\n", t)
		fmt.Printf("   - P4[anti-top] %v\n", tbar)
		fmt.Printf("\n")

		return nil
	})

	fmt.Printf("Number of events w/o reconstruction: %v\n\n", nBad)

	if err != nil {
		log.Fatalf("could not process tree: %+v", err)
	}

}
