# diagnostics_fill_state.F

First, add a block for the FFIELDS header under the GAD block (aronud line 28)
```
#ifdef ALLOW_ADDFLUID
#include "FFIELDS.h"
#endif
```

Second, add a line for addMass (around line 422):
```
#ifdef ALLOW_ADDFLUID
        CALL DIAGNOSTICS_FILL(addMass,'addMass ',0,Nr,0,1,1,myThid)
#endif /* ALLOW_ADDFLUID */
```
