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
		
		rvars = []rtree.ReadVar{
			{Name: "d_lep_pt" , Value: &lepPt},
			{Name: "d_lep_eta", Value: &lepEta},
			{Name: "d_lep_phi", Value: &lepPhi},
		}
	)

	// Get the TTree reader
	r, err := rtree.NewReader(t, rvars)
	if err != nil {
		log.Fatalf("could not create tree reader: %+v", err)
	}
	defer r.Close()

	// Event loop
	err = r.Read(func(ctx rtree.RCtx) error {
		fmt.Printf("evt[%d]: %v, %v, %v\n", ctx.Entry, lepPt, lepEta, lepPhi)
		return nil
	})
	if err != nil {
		log.Fatalf("could not process tree: %+v", err)
	}
	
}
