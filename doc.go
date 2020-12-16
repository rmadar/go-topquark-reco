// Package topreco allows to perform to top reconstruction in ttbar->dilepton events
// using analytical two methods. In a nutshell:
//
//  // Before the event loop
//  topBuilder, err := tbuilder.New("../testdata/smearingHistos.root")
//
//  // In the event loop: get the needed inputs
//  var (
//         // Lepton 4-momentum and flavour
//         lep0 fmom.PxPyPzE
//	   lep1 fmom.PxPyPzE
//         lid0 int
//	   lid1 int
//
//         // Jets 4-momentum and b-tagging info
//	   jet0 fmom.PxPyPzE
//	   jet1 fmom.PxPyPzE
//         j1b  bool
//	   j2b  bool
//
//         // Missing ET
//         etx  float64
//         ety  float64
//	)
//
//   // Call the reconstruction of top and antitop
//   t, tbar, nIterations := topBuilder.Reconstruct(
//	   lep0, lep1, lid0, lid1,
//	   jet0, jet1, j1b, j2b,
//	   etx, ety,
//      )
//
package topreco
