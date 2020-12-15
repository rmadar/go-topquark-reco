# How to use this tool in C/C++ with ROOT

### General methode

The user needs 3 files: `libtopreco.h`, `libtopreco.so` and `topreco.h`. In the ROOT macro, this tool can be called in the following way:

```cpp
#include "topreco.h"

// Initialization of the topBuilder (before the event loop)
auto smearFile = "../testdata/smearingHistos.root";
auto nSmearing = 100;
auto debugMode = 0;
InitTopBuilder((char*)smearFile, nSmearing, debugMode);

// Event loop
for (Int_t n = 0; n < N	; n++) {

    // Load the event
    t->GetEntry(n);

    // Get the relevant information for the reconstruction
    TLorentzVector l0, l1, j0, j1; // 4-momentum of two leptons and two jets
    int pid0, pid1                 // pid of lepton (used for a proper smearing)
    int j0b, j1b                   // btag info for the two jets
    double etx, ety                // missing transverse energy components

    // Define container to store the reconstrucred kinematics
    TLorentzVector t, tbar;

    // Run the reconstruction
    Int_t ok = runTopReconstruction(l0, l1, pid0, pid1,
       	       		            j0, j1, j0b, j1b,
				    etx, ety, &t, &tbar);
}
```

### Working example

A working example can be found [here](test-root.C)