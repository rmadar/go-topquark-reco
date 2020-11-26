package sonn

//#cgo LDFLAGS: -L/opt/root/root-6.22.02/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lROOTDataFrame -pthread -lm -ldl -rdynamic -lMathMore
//#cgo CXXFLAGS: -pthread -std=c++11 -m64 -I/opt/root/root-6.22.02/include
//
//#include "sonn.h"
//
//#include <stdlib.h> // for malloc/free
import "C"

import (
	"unsafe"

	"go-hep.org/x/hep/fmom"
)

func RecoTops(
	lep, lepbar fmom.PxPyPzE, pdgIDLep, pdgIDLepBar int,
	jet, jetbar fmom.PxPyPzE, emissx, emissy float64,
	nbjets int,
) []fmom.PxPyPzE {

	var (
		tlvLep = C.tlv_t{
			Px: C.double(lep.Px()),
			Py: C.double(lep.Py()),
			Pz: C.double(lep.Pz()),
			E:  C.double(lep.E()),
		}
		tlvLepBar = C.tlv_t{
			Px: C.double(lepbar.Px()),
			Py: C.double(lepbar.Py()),
			Pz: C.double(lepbar.Pz()),
			E:  C.double(lepbar.E()),
		}
		tlvJet = C.tlv_t{
			Px: C.double(jet.Px()),
			Py: C.double(jet.Py()),
			Pz: C.double(jet.Pz()),
			E:  C.double(jet.E()),
		}
		tlvJetBar = C.tlv_t{
			Px: C.double(jetbar.Px()),
			Py: C.double(jetbar.Py()),
			Pz: C.double(jetbar.Pz()),
			E:  C.double(jetbar.E()),
		}
	)

	out := C.sonn(
		tlvLep, tlvLepBar, C.int(pdgIDLep), C.int(pdgIDLepBar),
		tlvJet, tlvJetBar, C.double(emissx), C.double(emissy),
		C.int(nbjets),
	)
	defer C.free(unsafe.Pointer(out.tlvs))

	tlvs := make([]fmom.PxPyPzE, out.n)
	for i := range tlvs {
		tlv := C.get_tlv(out.tlvs, C.int(i))
		tlvs[i] = fmom.NewPxPyPzE(
			float64(tlv.Px),
			float64(tlv.Py),
			float64(tlv.Pz),
			float64(tlv.E),
		)
	}

	return tlvs
}

/* TO BE WRITTEN
func SetupSmearingFile(fname string) {
	C.sonn()
}
*/
