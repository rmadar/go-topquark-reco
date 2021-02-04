# Top Quark Reconstruction

This repository holds a pure go implementation of the Sonneschein method ([arXiv:hep-ph/0603011](https://arxiv.org/abs/hep-ph/0603011))
to reconstruct top quark in di-lepton events. The documentation of the package can be found [here](https://godoc.org/github.com/rmadar/go-topquark-reco/sonn). A wraper to use this tool in C is under developpement, enabeling its use in a ROOT macro (cf. [capi](capi)).

## Workflow of the method/code

This tool proceeds with the following steps for each event:
 1. Loop over the jets combinatorics (currently, `(0, 1)` and `(1, 0)` are considered)
 2. Produce `Nsmear` smeared kinematics based on the PDFs of relevant variables (stored in `smearingHistos.root`)
 3. Call the recontruction(s) for each of these smeared events, average the result and assign a final weight
 4. Compare which jet configuration give the highest weight and select it.

## In a nutshell

Before the event loop starts, the tool must be initialized (few options are possibles) with the line

```go
 // Before the event loop
 topBuilder, err := tbuilder.New("../testdata/smearingHistos.root")
```

For each event, the relevant four-momentum must be gathered and three
methods can be called depending on which methods you want. **Important note:** the third
method `AllReco()` is faster than calling both `SonnReco()` and `ElliReco()`
because the smearing of the kinematics if performed only once (while it would be done twice if
the two individual functions would be called). In that case, the returned objects
are four-momentum arrays of 2 elements each (`ts[0]` is the top 4-momentum obtained by the
Sonnenschein methods, while `ts[1]` is the top 4-momentum obtained by the Ellipse method).
 
 
```go
 // In the event loop: get the needed inputs
 var (
        // Lepton 4-momentum and flavour
        lep0, lep1 fmom.PxPyPzE
        lid0, lid1 int

        // Jets 4-momentum and b-tagging info
        jet0, jet1 fmom.PxPyPzE
        j1b , j2b  bool

        // Missing ET
        etx, ety  float64
     )

  // Call the Sonnenschein reconstruction only
  t, tbar, nIterations := topBuilder.SonnReco(
	   lep0, lep1, lid0, lid1,
	   jet0, jet1, j1b, j2b,
	   etx, ety,
     )

  // Call the Ellipse reconstruction only
  t, tbar, nIterations := topBuilder.ElliReco(
	   lep0, lep1, lid0, lid1,
	   jet0, jet1, j1b, j2b,
	   etx, ety,
     )

  // Call the both reconstructions, avoiding to perform the smearing twice.
  ts, tbars, nIterations := topBuilder.AllReco(
	   lep0, lep1, lid0, lid1,
	   jet0, jet1, j1b, j2b, 
	   etx, ety,   
     )	   	

```