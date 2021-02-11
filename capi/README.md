# How to use this tool in C/C++ with ROOT

**THIS SECTION IS CURRENTLY OUT-OF-DATE**

### General methode

The user needs 3 files: `libtopreco.h`, `libtopreco.so` and `topreco.h`. In the ROOT macro, this tool can be called in the following way:

```cpp
#include "topreco.h"

// Initialization the Top Reconstruction tool (before the event loop)
auto smearFile = "../testdata/smearingHistos.root";
auto nSmearing = 100;
auto debugMode = false;
initTopReconstruction(smearFile, nSmearing, debugMode);

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

    // Run the reconstruction, return the number of successful
    // smearing iterations averaged to produce the final kinematics.
    Int_t nIter = runTopReconstruction(l0, l1, pid0, pid1,
       	       		               j0, j1, j0b, j1b,
				       etx, ety, &t, &tbar);
}
```

### Working example

A working example can be found [here](test-root.C)