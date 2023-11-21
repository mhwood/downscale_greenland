# apply_forcing.F functions

The changes in the apply_forcing.F script are in two functions: APPLY_FORCING_T and APPLY_FORCING_S

## Changes to APPLY_FORCING_T

Inside the APPLY_FORCING_T, the iceplume header files is included (around line 415):
```
#ifdef ALLOW_ICEPLUME
#include "ICEPLUME.h"
#endif /* ALLOW_ICEPLUME */
```
The tmpVar is zero'd after the FIZHI tendency apply (aroudn line 594)
```
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
        tmpVar(i,j)=0. _d 0
       ENDDO
      ENDDO
```
Then, inside and the start of the #ifdef ALLOW_ADDFLUID block, the tendencies are applied to gT_arr:
```
#ifdef ALLOW_ICEPLUME

      IF ( useICEPLUME ) THEN
       IF ( selectAddFluid.NE.0) THEN
        IF ( ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &      .OR. convertFW2Salt.EQ.-1. _d 0 ) THEN
         DO j=0,sNy+1
          DO i=0,sNx+1
           tmpVar(i,j)=
     &          addMass3Dplume(i,j,k,bi,bj)*mass2rUnit
     &          *(temp_addMass3Dplume(I,J,k,bi,bj)-theta(i,j,k,bi,bj))
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)

           gT_arr(i,j) = gT_arr(i,j) + tmpVar(i,j)
          ENDDO
         ENDDO
       ELSE
         DO j=0,sNy+1
          DO i=0,sNx+1
           tmpVar(i,j)=
     &          addMass3Dplume(i,j,k,bi,bj)*mass2rUnit
     &          *( temp_addMass3Dplume(I,J,k,bi,bj) - tRef(k) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)

           gT_arr(i,j) = gT_arr(i,j) + tmpVar(i,j)
          ENDDO
         ENDDO
        ENDIF
       ENDIF
      ENDIF
#else /* ALLOW_ICEPLUME */
```
Then, the existing loops inside the ALLOW_ADDFLUID block are updated:
```
      IF ( selectAddFluid.NE.0 .AND. temp_addMass.NE.UNSET_RL ) THEN
       IF ( ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &      .OR. convertFW2Salt.EQ.-1. _d 0 ) THEN
         DO j=0,sNy+1
          DO i=0,sNx+1

            tmpVar(i,j) =
     &          addMass(i,j,k,bi,bj)*mass2rUnit
     &          *( temp_addMass - theta(i,j,k,bi,bj) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)

            gT_arr(i,j) = gT_arr(i,j) + tmpVar(i,j)
          ENDDO
         ENDDO
       ELSE
         DO j=0,sNy+1
          DO i=0,sNx+1
            tmpVar(i,j) = 
     &          addMass(i,j,k,bi,bj)*mass2rUnit
     &          *( temp_addMass - tRef(k) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)

            gT_arr(i,j) = gT_arr(i,j) + tmpVar(i,j)
          ENDDO
         ENDDO
       ENDIF
      ENDIF
```
Finally, still inside the ALLOWLADDFLUID block, a new diagnostics is added:
```
#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
          CALL DIAGNOSTICS_FILL(tmpVar,
     &                     'gTaddMas',k,1,2,bi,bj,myThid )
#ifdef ALLOW_ICEPLUME
        DO j=0,sNy+1
          DO i=0,sNx+1
            tmpVar(i,j)=theta(i,j,k,bi,bj)
     &          * addMass3Dplume(i,j,k,bi,bj)*mass2rUnit
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
          ENDDO
        ENDDO
          CALL DIAGNOSTICS_FILL(tmpVar,
     &                     'tAddMass',k,1,2,bi,bj,myThid )
#endif /* ALLOW_ICEPLUME */
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */
```
After all of the tendencies have been updates, the tendency is applied after the icefront tendency apply (arounf line 836):
```
#ifdef ALLOW_ICEPLUME
      IF ( useICEPLUME )
     &     CALL ICEPLUME_TENDENCY_APPLY_T(
     U                   gT_arr,
     I                   iMin,iMax,jMin,jMax, 
     I                   k, bi, bj, myTime, myIter, myThid )
#endif /* ALLOW_ICEPLUME */
```

## Changes to APPLY_FORCING_S

Inside the APPLY_FORCING_S, the iceplume header files is included (around line 916):
```
#ifdef ALLOW_ICEPLUME
#include "ICEPLUME.h"
#endif /* ALLOW_ICEPLUME */
```
Then, a tmpVar is defined in the header block (around line 938):
```
      _RL     tmpVar(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
```
Then, just as for temperature, the tmpVar is zero'd after the FIZHI_TENDENCY_APPLY function:
```
      DO j=1-OLy,sNy+OLy
       DO i=1-OLx,sNx+OLx
        tmpVar(i,j) = 0. _d 0
       ENDDO
      ENDDO
```
Similarly, the ALLOW_ADDFLUID block is updated as follows:
```
      IF ( useICEPLUME ) THEN
       IF ( selectAddFluid.NE.0) THEN
        IF ( ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &      .OR. convertFW2Salt.EQ.-1. _d 0 ) THEN
         DO j=0,sNy+1
          DO i=0,sNx+1

           tmpVar(i,j) =
     &          addMass3Dplume(i,j,k,bi,bj)*mass2rUnit
     &          *( salt_addMass3Dplume(I,J,k,bi,bj)-salt(i,j,k,bi,bj) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)

           gS_arr(i,j) = gS_arr(i,j) + tmpVar(i,j)
          ENDDO
         ENDDO
       ELSE
         DO j=0,sNy+1
          DO i=0,sNx+1

           tmpVar(i,j) =
     &          addMass3Dplume(i,j,k,bi,bj)*mass2rUnit
     &          *( salt_addMass3Dplume(I,J,k,bi,bj) - sRef(k) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)

           gS_arr(i,j) = gS_arr(i,j) + tmpVar(i,j)
          ENDDO
         ENDDO
        ENDIF
       ENDIF
      ENDIF
#else
catn: make sure (salt,temp)_addMass are unset so we dont enter here
      IF ( selectAddFluid.NE.0 .AND. salt_addMass.NE.UNSET_RL ) THEN
       IF ( ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &      .OR. convertFW2Salt.EQ.-1. _d 0 ) THEN

         DO j=0,sNy+1
          DO i=0,sNx+1

            tmpVar(i,j) =
     &          addMass(i,j,k,bi,bj)*mass2rUnit
     &          *( salt_addMass - salt(i,j,k,bi,bj) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)

            gS_arr(i,j) = gS_arr(i,j) + tmpVar(i,j)
          ENDDO
         ENDDO
       ELSE
         DO j=0,sNy+1
          DO i=0,sNx+1

            tmpVar(i,j) =
     &          addMass(i,j,k,bi,bj)*mass2rUnit
     &          *( salt_addMass - sRef(k) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)

            gS_arr(i,j) = gS_arr(i,j) + tmpVar(i,j)
          ENDDO
         ENDDO
       ENDIF
      ENDIF
#endif /* ALLOW_ICEPLUME */

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
          CALL DIAGNOSTICS_FILL(tmpVar,
     &                     'gSaddMas',k,1,2,bi,bj,myThid )
#ifdef ALLOW_ICEPLUME
         DO j=0,sNy+1
          DO i=0,sNx+1
            tmpVar(i,j)=salt(i,j,k,bi,bj)
     &          *addMass3Dplume(i,j,k,bi,bj)*mass2rUnit
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(k)*_recip_hFacC(i,j,k,bi,bj)
          ENDDO
        ENDDO
          CALL DIAGNOSTICS_FILL(tmpVar,
     &                     'sAddMass',k,1,2,bi,bj,myThid )
#endif /* ALLOW_ICEPLUME */
      ENDIF
#endif /* ALLOW_DIAGNOSTICS */
```
The last change is to call the TENDENCY_APPLY function after the same one for icefront:
```
#ifdef ALLOW_ICEPLUME
      IF ( useICEPLUME )
     &     CALL ICEPLUME_TENDENCY_APPLY_S(
     U                   gS_arr,
     I                   iMin,iMax,jMin,jMax,
     I                   k, bi, bj, myTime, myIter, myThid )
#endif /* ALLOW_ICEPLUME */
```
