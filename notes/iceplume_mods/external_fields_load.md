# external_fields_load.md

This script is modified with the following components:

1. A header include block added under the `#include "DYNVARS.h"` block:
```
#ifdef ALLOW_ICEPLUME
#include "ICEPLUME.h"
#endif /* ALLOW_ICEPLUME */
```

2. An ALLOW_ICEPLUME block after the `#ifdef ATMOSPHERIC_LOADING` block around line 211:
```
# ifdef ALLOW_ICEPLUME
      IF ( useICEPLUME ) THEN
      IF ( runoffQsgfile .NE. ' ') THEN
      CALL READ_REC_XY_RL(
     &     runoffQsgFile,runoffQsg0,inTime0,myIter,myThid)
      CALL READ_REC_XY_RL(
     &     runoffQsgFile,runoffQsg1,inTime1,myIter,myThid)
      ENDIF
      ENDIF
# endif /* ALLOW_ICEPLUME */
```

3. Another ALLOW_ICEPLUME block after the `ATMOSPHERIC_LOADING` block around line 350:
```
#ifdef ALLOW_ICEPLUME
        IF ( runoffQsgFile .NE. ' '  ) THEN
          DO j=1-Oly,sNy+Oly
           DO i=1-Olx,sNx+Olx
            runoffQsg(i,j,bi,bj) = bWght*runoffQsg0(i,j,bi,bj)
     &                        + aWght*runoffQsg1(i,j,bi,bj)
           ENDDO
          ENDDO
        ENDIF
#endif
```
