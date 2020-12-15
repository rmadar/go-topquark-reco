# How to use this tool in C/C++ with ROOT

### General methode

The user needs 3 files: `libtopreco.h`, `libtopreco.so` and `topreco.h`. In the ROOT macro, this tool can be called in the following way:

```cpp
#include "topreco.h"

// Initialization of the topBuilder (before the event loop)
InitTopBuilder((char*)"../testdata/smearingHistos.root");

// Event loop
for (Int_t n = 0; n < N	; n++) {

    // Load the event
    t->GetEntry(n);

    // Get the 4-momentum, pid, b-tag and missing ET
    // ...

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