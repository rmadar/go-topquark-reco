# Top Quark Reconstruction

This repository holds a pure go implementation of the Sonneschein method ([arXiv:hep-ph/0603011](https://arxiv.org/abs/hep-ph/0603011))
to reconstruct top quark in di-lepton events. The documentation of the package can be found [here](https://godoc.org/github.com/rmadar/go-topquark-reco/sonn). A wraper to use this tool in C is under developpement, enabeling its use in a ROOT macro (cf. [capi](capi)).

## In a nutshell

```go
 // Before the event loop
 topBuilder, err := tbuilder.New("../testdata/smearingHistos.root")

 // In the event loop: get the needed inputs
 var (
        // Lepton 4-momentum and flavour
        lep0 fmom.PxPyPzE
        lep1 fmom.PxPyPzE
        lid0 int
        lid1 int

        // Jets 4-momentum and b-tagging info
        jet0 fmom.PxPyPzE
        jet1 fmom.PxPyPzE
        j1b  bool
        j2b  bool

        // Missing ET
        etx  float64
        ety  float64
     )

  // Call the Sonnenschein reconstruction
  t, tbar, nIterations := topBuilder.SonnReco(
	   lep0, lep1, lid0, lid1,
	   jet0, jet1, j1b, j2b,
	   etx, ety,
     )

  // Call the Ellipse reconstruction
  t, tbar, nIterations := topBuilder.ElliReco(
	   lep0, lep1, lid0, lid1,
	   jet0, jet1, j1b, j2b,
	   etx, ety,
     )

  // Call the both reconstructions - avoiding to perform the smearing twice.
  ts, tbars, nIterations := topBuilder.ElliReco(
	   lep0, lep1, lid0, lid1,
	   jet0, jet1, j1b, j2b, 
	   etx, ety,   
     )	   	

```