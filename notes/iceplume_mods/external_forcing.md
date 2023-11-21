# external_forcing.md

The changes to this script are made in the EXTERNAL_FORCING_T and EXTERNAL_FORCING_S functions:

## EXTERNAL_FORCING_T
Add an include block below the SURFACE.h include (around line 316):
```
#ifdef ALLOW_ICEPLUME
#include "ICEPLUME.h"
#endif /* ALLOW_ICEPLUME */
```

Update the ALLOW_ADDFLUID (around line 380) to be:
```
#ifdef ALLOW_ADDFLUID

#ifdef ALLOW_ICEPLUME
      DO j=1,sNy
       DO i=1,sNx
        tmpVar(i,j)=0. _d 0
       ENDDO
      ENDDO

      IF ( useICEPLUME ) THEN
       IF ( selectAddFluid.NE.0 ) THEN
        IF ( ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &      .OR. convertFW2Salt.EQ.-1. _d 0 ) THEN
         DO j=1,sNy
          DO i=1,sNx
            tmpVar(i,j) =
     &        addMass3Dplume(i,j,kLev,bi,bj)*mass2rUnit
     &        *(temp_addMass3D(I,J,Klev,bi,bj)-theta(i,j,kLev,bi,bj) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)

            gT(i,j,klev,bi,bj) = gT(i,j,klev,bi,bj) + tmpVar(i,j)

          ENDDO
         ENDDO
        ELSE
         DO j=1,sNy
          DO i=1,sNx
            tmpVar(i,j) = 
     &        addMass3Dplume(i,j,kLev,bi,bj)*mass2rUnit
     &          *( temp_addMass3D(I,J,Klev,bi,bj) - tRef(kLev) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)

            gT(i,j,klev,bi,bj) = gT(i,j,klev,bi,bj) + tmpVar(i,j)
          ENDDO
         ENDDO
        ENDIF
       ENDIF
      ENDIF

#else
      IF ( selectAddFluid.NE.0 .AND. temp_addMass.NE.UNSET_RL ) THEN
       IF ( ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &      .OR. convertFW2Salt.EQ.-1. _d 0 ) THEN
         DO j=1,sNy
          DO i=1,sNx

            tmpVar(i,j) =
     &          addMass(i,j,kLev,bi,bj)*mass2rUnit
     &          *( temp_addMass - theta(i,j,kLev,bi,bj) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)

            gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj) + tmpVar(i,j)
          ENDDO
         ENDDO
       ELSE
         DO j=1,sNy
          DO i=1,sNx

            tmpVar(i,j) =
     &          addMass(i,j,kLev,bi,bj)*mass2rUnit
     &          *( temp_addMass - tRef(kLev) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)

            gT(i,j,kLev,bi,bj) = gT(i,j,kLev,bi,bj) + tmpVar(i,j)
          ENDDO
         ENDDO
       ENDIF
      ENDIF
#endif /* ALLOW_ICEPLUME */
#endif /* ALLOW_ADDFLUID */
```
Then, call the ICEPLUME_TENDENCY_APPLY after the ALLOW_ICEFRONT block (around line 573):
```
#ifdef ALLOW_ICEPLUME
      IF ( useICEPLUME )
     &     CALL ICEPLUME_TENDENCY_APPLY_T(
     U                   gT(1-OLx,1-OLy,kLev,bi,bj),
     I                   iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                   kLev, bi, bj, myTime, 0, myThid )
#endif /* ALLOW_ICEPLUME */
```


## EXTERNAL_FORCING_S
Add an include block below the SURFACE.h include (around line 714):
```
#ifdef ALLOW_ICEPLUME
#include "ICEPLUME.h"
#endif /* ALLOW_ICEPLUME */
```

Update the ALLOW_ADDFLUID (around line 702) to be:
```
#ifdef ALLOW_ADDFLUID

#ifdef ALLOW_ICEPLUME
      DO j=0,sNy+1
       DO i=0,sNx+1
        tmpVar(i,j)=0. _d 0
       ENDDO
      ENDDO

      IF ( useICEPLUME ) THEN
       IF ( selectAddFluid.NE.0 ) THEN
        IF ( ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &      .OR. convertFW2Salt.EQ.-1. _d 0 ) THEN
         DO j=1,sNy
          DO i=1,sNx
            tmpVar(i,j) = 
     &        addMass3Dplume(i,j,kLev,bi,bj)*mass2rUnit
     &        *( salt_addMass3D(I,J,Klev,bi,bj) - salt(i,j,kLev,bi,bj) )
     &        *recip_rA(i,j,bi,bj)
     &        *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
c     &        *recip_deepFac2C(kLev)*recip_rhoFacC(kLev)

            gS(i,j,klev,bi,bj)=gS(i,j,k,bi,bj) + tmpVar(i,j)

          ENDDO
         ENDDO
        ELSE
         DO j=1,sNy
          DO i=1,sNx

            tmpVar(i,j) =
     &        addMass3Dplume(i,j,kLev,bi,bj)*mass2rUnit
     &          *( salt_addMass3D(I,J,Klev,bi,bj) - sRef(kLev) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
c     &         *recip_deepFac2C(kLev)*recip_rhoFacC(kLev)

            gS(i,j,klev,bi,bj)=gS(i,j,klev,bi,bj) + tmpVar(i,j)

          ENDDO
         ENDDO
        ENDIF
       ENDIF
      ENDIF
#else
      IF ( selectAddFluid.NE.0 .AND. salt_addMass.NE.UNSET_RL ) THEN
       IF ( ( selectAddFluid.GE.1 .AND. nonlinFreeSurf.GT.0 )
     &      .OR. convertFW2Salt.EQ.-1. _d 0 ) THEN
         DO j=1,sNy
          DO i=1,sNx

            tmpVar(i,j) =
     &          addMass(i,j,kLev,bi,bj)*mass2rUnit
     &          *( salt_addMass - salt(i,j,kLev,bi,bj) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
C    &          *recip_deepFac2C(kLev)*recip_rhoFacC(kLev)

            gS(i,j,kLev,bi,bj) = gS(i,j,kLev,bi,bj) + tmpVar(i,j)
          ENDDO
         ENDDO
       ELSE
         DO j=1,sNy
          DO i=1,sNx

            tmpVar(i,j) =
     &          addMass(i,j,kLev,bi,bj)*mass2rUnit
     &          *( salt_addMass - sRef(kLev) )
     &          *recip_rA(i,j,bi,bj)
     &          *recip_drF(kLev)*_recip_hFacC(i,j,kLev,bi,bj)
C    &          *recip_deepFac2C(kLev)*recip_rhoFacC(kLev)

            gS(i,j,kLev,bi,bj) = gS(i,j,kLev,bi,bj) + tmpVar(i,j)
          ENDDO
         ENDDO
       ENDIF
      ENDIF
#endif /* ALLOW_ICEPLUME */
#endif /* ALLOW_ADDFLUID */
```
Then, call the ICEPLUME_TENDENCY_APPLY after the ALLOW_ICEFRONT block (around line 901):
```
#ifdef ALLOW_ICEPLUME
      IF ( useICEPLUME )
     &     CALL ICEPLUME_TENDENCY_APPLY_S(
     U                   gS(1-OLx,1-OLy,kLev,bi,bj),
     I                   iMin,iMax,jMin,jMax, kLev, bi,bj,
     I                   kLev, bi, bj, myTime, 0, myThid )
#endif /* ALLOW_ICEPLUME */
```

