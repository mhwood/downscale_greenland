# exf_getffields.F

In this script, an include statement is made for the ICEPLUME header under the `#include "EXF_FIELDS.h"` line (around line 41):
```
#ifdef ALLOW_ICEPLUME
#include "ICEPLUME.h"
#endif /* ALLOW_ICEPLUME */
```
And a block for ICEPLUME is added after the ALLOW_RUNOFFTEMP block:
```
#ifdef ALLOW_ICEPLUME
      IF (useICEPLUME) THEN
      CALL EXF_SET_FLD(
     I     'runoffQsg', runoffQsgfile, runoffQsgmask,
     I     runoffQsgStartTime, runoffQsgperiod, runoffQsgRepCycle,
     I     runoffQsg_inscal,
     I     runoffQsg_remov_intercept, runoffQsg_remov_slope,
     U     runoffQsg, runoffQsg0, runoffQsg1,
#ifdef USE_EXF_INTERPOLATION
     I     runoffQsg_lon0, runoffQsg_lon_inc, 
     I     runoffQsg_lat0, runoffQsg_lat_inc,
     I     runoffQsg_nlon, runoffQsg_nlat, xC, yC, 
     I     runoffQsg_interpMethod,
#endif
     I     myTime, myIter, myThid )
      ENDIF
#endif /* ALLOW_ICEPLUME */
```
