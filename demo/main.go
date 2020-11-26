package main

import (
	"log"
	"fmt"

	"go-hep.org/x/hep/groot"
	"go-hep.org/x/hep/groot/rtree"
	
	_ "github.com/rmadar/go-topquark-reco/sonn"
)

func main(){

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
		lepPt  []float32
		lepEta []float32
		lepPhi []float32
		lepPid []int32
		jetPt  []float32
		jetEta []float32
		jetPhi []float32
		nBjets   int32
		metMet   float32
		metPhi   float32
		
		rvars = []rtree.ReadVar{
			{Name: "d_lep_pt" , Value: &lepPt},
			{Name: "d_lep_eta", Value: &lepEta},
			{Name: "d_lep_phi", Value: &lepPhi},
			{Name: "d_lep_pid", Value: &lepPid},
			{Name: "d_jet_pt" , Value: &jetPt},
			{Name: "d_jet_eta", Value: &jetEta},
			{Name: "d_jet_phi", Value: &jetPhi},
			{Name: "d_nbjet"  , Value: &nBjets},
			{Name: "d_met_met", Value: &metMet},
			{Name: "d_met_phi", Value: &metPhi},
		}
	)

	// Get the TTree reader
	r, err := rtree.NewReader(t, rvars)
	if err != nil {
		log.Fatalf("could not create tree reader: %+v", err)
	}
	defer r.Close()

	// Load smearing histograms
	// sonn.SetupSmearingFile("")
	
	// Event loop
	err = r.Read(func(ctx rtree.RCtx) error {

		// Print some information
		fmt.Printf("evt[%d]: %v, %v, %v\n", ctx.Entry, lepPt, lepEta, lepPhi)

		// Call the Sonnenschein reconstruction
		// translate input TTree variable into function argument
		//   * lP, lbarP, jP, jbarP -> fmom.PxPyPzE()
		//   * lId, lbId, nBjets    -> int
		//   * Etx, Ety             -> float32 
		//tops := sonn.RecoTops(lP, lbarP, lId, lbarId, jP, jbarP, Etx, Ety, nBjets)
		
		return nil
	})
	if err != nil {
		log.Fatalf("could not process tree: %+v", err)
	}
	
}
