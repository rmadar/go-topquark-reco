package main

// typedef struct {
//     double Px, Py, Pz, E;
// } P4_t;
import "C"
import (
	"log"

	"github.com/rmadar/go-topquark-reco/tbuilder"
	"go-hep.org/x/hep/fmom"
)

func main() {}

var bldr *tbuilder.TopBuilder

//export InitTopBuilder
func InitTopBuilder(fname *C.char, nSmear C.int, debug C.int) {
	
	debugB := false
	if debug == 1 {
		debugB = true
	}

	var err error
	bldr, err = tbuilder.New(C.GoString(fname),
		tbuilder.WithSmearN(int(nSmear)),
		tbuilder.WithDebug(debugB),
	)
	
	if err != nil {
		log.Panicf("could not create TopBuilder: %+v", err)
	}
}

//export ReconstructTops
func ReconstructTops(
	lep0, lep1 C.P4_t, lep0PDG, lep1PDG int,
	jet0, jet1 C.P4_t, bjet0, bjet1 int,
	emissx, emissy float64,
	top0, top1 *C.P4_t) int {

	*top0 = C.P4_t{}
	*top1 = C.P4_t{}

	t0, t1, ok := bldr.Reconstruct(
		newPxPyPzE(lep0), newPxPyPzE(lep1), lep0PDG, lep1PDG,
		newPxPyPzE(jet0), newPxPyPzE(jet1), bjet0 != 0, bjet1 != 0,
		float64(emissx), float64(emissy),
	)

	if ok != 0 {
		*top0 = C.P4_t{C.double(t0.Px()), C.double(t0.Py()), C.double(t0.Pz()), C.double(t0.E())}
		*top1 = C.P4_t{C.double(t1.Px()), C.double(t1.Py()), C.double(t1.Pz()), C.double(t1.E())}
	}
	return ok
}

func newPxPyPzE(p4 C.P4_t) fmom.PxPyPzE {
	return fmom.NewPxPyPzE(
		float64(p4.Px), float64(p4.Py), float64(p4.Pz),
		float64(p4.E),
	)
}

func newP4(p4 fmom.PxPyPzE) C.P4_t {
	return C.P4_t{
		C.double(p4.Px()),
		C.double(p4.Py()),
		C.double(p4.Pz()),
		C.double(p4.E()),
	}
}
