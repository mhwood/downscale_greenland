#include "DARWIN_OPTIONS.h"

CBOP
C !ROUTINE: DARWIN_SINKING
C !INTERFACE: ==========================================================
      SUBROUTINE DARWIN_SINKING(
     I     Ptr,
     U     gTr,
     I     bi,bj,iMin,iMax,jMin,jMax,myTime,myIter,myThid)

C !DESCRIPTION:
C     compute tendencies from sinking of particulate organic matter

C !USES: ===============================================================
      IMPLICIT NONE
#include "SIZE.h"
c#include "EEPARAMS.h"
#include "GRID.h"
#include "DARWIN_SIZE.h"
#include "DARWIN_INDICES.h"
#include "DARWIN_PARAMS.h"
#include "DARWIN_TRAITS.h"

C !INPUT PARAMETERS: ===================================================
C  myThid :: thread number
C  Ptr    :: darwin model tracers
C  dT     :: time step (for fake tendency with useEquilibriumC)
      _RL Ptr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy,nDarwin)
      INTEGER bi,bj,iMin,iMax,jMin,jMax
      INTEGER myThid, myIter
      _RL myTime

C !INPUT/OUTPUT PARAMETERS: ============================================
C  gTr    :: computed tendencies
      _RL gTr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nDarwin)
CEOP

#ifdef ALLOW_DARWIN

c !LOCAL VARIABLES: ====================================================
      INTEGER i,j,k,l
      _RL dzup, dzdn, flux

      DO k=1,Nr-1
       DO j=jMin,jMax
        DO i=iMin,iMax
         dzup = DRF(k)*hFacC(i,j,k,bi,bj)
         dzdn = DRF(k+1)*hFacC(i,j,k+1,bi,bj)
         IF (dzdn .GT. 0 _d 0 .AND. dzup .GT. 0 _d 0) THEN
          flux = wPIC_sink*MAX(0.0, Ptr(i,j,k,bi,bj,iPIC))
          gTr(i,j,k  ,iPIC ) = gTr(i,j,k  ,iPIC ) - flux/dzup
          gTr(i,j,k+1,iPIC ) = gTr(i,j,k+1,iPIC ) + flux/dzdn
          flux = wC_sink*MAX(0.0, Ptr(i,j,k,bi,bj,iPOC))
          gTr(i,j,k  ,iPOC ) = gTr(i,j,k  ,iPOC ) - flux/dzup
          gTr(i,j,k+1,iPOC ) = gTr(i,j,k+1,iPOC ) + flux/dzdn
          flux = wN_sink*MAX(0.0, Ptr(i,j,k,bi,bj,iPON))
          gTr(i,j,k  ,iPON ) = gTr(i,j,k  ,iPON ) - flux/dzup
          gTr(i,j,k+1,iPON ) = gTr(i,j,k+1,iPON ) + flux/dzdn
          flux = wP_sink*MAX(0.0, Ptr(i,j,k,bi,bj,iPOP))
          gTr(i,j,k  ,iPOP ) = gTr(i,j,k  ,iPOP ) - flux/dzup
          gTr(i,j,k+1,iPOP ) = gTr(i,j,k+1,iPOP ) + flux/dzdn
          flux = wSi_sink*MAX(0.0, Ptr(i,j,k,bi,bj,iPOSi))
          gTr(i,j,k  ,iPOSi) = gTr(i,j,k  ,iPOSi) - flux/dzup
          gTr(i,j,k+1,iPOSi) = gTr(i,j,k+1,iPOSi) + flux/dzdn
          flux = wFe_sink*MAX(0.0, Ptr(i,j,k,bi,bj,iPOFe))
          gTr(i,j,k  ,iPOFe) = gTr(i,j,k  ,iPOFe) - flux/dzup
          gTr(i,j,k+1,iPOFe) = gTr(i,j,k+1,iPOFe) + flux/dzdn
          DO l = 1, nplank
           flux = biosink(l)*MAX(0.0, Ptr(i,j,k,bi,bj,ic+l-1))
           gTr(i,j,k  ,ic+l-1 )=gTr(i,j,k  ,ic+l-1 ) - flux/dzup
           gTr(i,j,k+1,ic+l-1 )=gTr(i,j,k+1,ic+l-1 ) + flux/dzdn
           flux = bioswim(l)*MAX(0.0, Ptr(i,j,k+1,bi,bj,ic+l-1))
           gTr(i,j,k  ,ic+l-1 )=gTr(i,j,k  ,ic+l-1 ) + flux/dzup
           gTr(i,j,k+1,ic+l-1 )=gTr(i,j,k+1,ic+l-1 ) - flux/dzdn
#ifdef DARWIN_ALLOW_NQUOTA
           flux = biosink(l)*MAX(0.0, Ptr(i,j,k,bi,bj,in+l-1))
           gTr(i,j,k  ,in+l-1 )=gTr(i,j,k  ,in+l-1 ) - flux/dzup
           gTr(i,j,k+1,in+l-1 )=gTr(i,j,k+1,in+l-1 ) + flux/dzdn
           flux = bioswim(l)*MAX(0.0, Ptr(i,j,k+1,bi,bj,in+l-1))
           gTr(i,j,k  ,in+l-1 )=gTr(i,j,k  ,in+l-1 ) + flux/dzup
           gTr(i,j,k+1,in+l-1 )=gTr(i,j,k+1,in+l-1 ) - flux/dzdn
#endif
#ifdef DARWIN_ALLOW_PQUOTA
           flux = biosink(l)*MAX(0.0, Ptr(i,j,k,bi,bj,ip+l-1))
           gTr(i,j,k  ,ip+l-1 )=gTr(i,j,k  ,ip+l-1 ) - flux/dzup
           gTr(i,j,k+1,ip+l-1 )=gTr(i,j,k+1,ip+l-1 ) + flux/dzdn
           flux = bioswim(l)*MAX(0.0, Ptr(i,j,k+1,bi,bj,ip+l-1))
           gTr(i,j,k  ,ip+l-1 )=gTr(i,j,k  ,ip+l-1 ) + flux/dzup
           gTr(i,j,k+1,ip+l-1 )=gTr(i,j,k+1,ip+l-1 ) - flux/dzdn
#endif
#ifdef DARWIN_ALLOW_SIQUOTA
           flux = biosink(l)*MAX(0.0, Ptr(i,j,k,bi,bj,isi+l-1))
           gTr(i,j,k  ,isi+l-1)=gTr(i,j,k  ,isi+l-1) - flux/dzup
           gTr(i,j,k+1,isi+l-1)=gTr(i,j,k+1,isi+l-1) + flux/dzdn
           flux = bioswim(l)*MAX(0.0, Ptr(i,j,k+1,bi,bj,isi+l-1))
           gTr(i,j,k  ,isi+l-1)=gTr(i,j,k  ,isi+l-1) + flux/dzup
           gTr(i,j,k+1,isi+l-1)=gTr(i,j,k+1,isi+l-1) - flux/dzdn
#endif
#ifdef DARWIN_ALLOW_FEQUOTA
           flux = biosink(l)*MAX(0.0, Ptr(i,j,k,bi,bj,ife+l-1))
           gTr(i,j,k  ,ife+l-1)=gTr(i,j,k  ,ife+l-1) - flux/dzup
           gTr(i,j,k+1,ife+l-1)=gTr(i,j,k+1,ife+l-1) + flux/dzdn
           flux = bioswim(l)*MAX(0.0, Ptr(i,j,k+1,bi,bj,ife+l-1))
           gTr(i,j,k  ,ife+l-1)=gTr(i,j,k  ,ife+l-1) + flux/dzup
           gTr(i,j,k+1,ife+l-1)=gTr(i,j,k+1,ife+l-1) - flux/dzdn
#endif
          ENDDO
#ifdef DARWIN_ALLOW_CHLQUOTA
          DO l = 1, nPhoto
           flux = biosink(l)*MAX(0.0, Ptr(i,j,k,bi,bj,iChl+l-1))
           gTr(i,j,k  ,iChl+l-1)=gTr(i,j,k  ,iChl+l-1)-flux/dzup
           gTr(i,j,k+1,iChl+l-1)=gTr(i,j,k+1,iChl+l-1)+flux/dzdn
           flux = bioswim(l)*MAX(0.0, Ptr(i,j,k+1,bi,bj,iChl+l-1))
           gTr(i,j,k  ,iChl+l-1)=gTr(i,j,k  ,iChl+l-1)+flux/dzup
           gTr(i,j,k+1,iChl+l-1)=gTr(i,j,k+1,iChl+l-1)-flux/dzdn
          ENDDO
#endif
         ENDIF
        ENDDO
       ENDDO
      ENDDO

#ifdef DARWIN_BOTTOM_SINK
      DO k=1,Nr
       DO j=jMin,jMax
        DO i=iMin,iMax
         dzup = DRF(k)*hFacC(i,j,k,bi,bj)
         IF (k.LT.Nr) THEN
           dzdn = DRF(k+1)*hFacC(i,j,k+1,bi,bj)
         ELSE
           dzdn = 0. _d 0
         ENDIF
         IF (dzup .GT. 0 _d 0 .AND. dzdn .EQ. 0 _d 0) THEN
          flux = wPIC_sink*MAX(0.0, Ptr(i,j,k,bi,bj,iPIC))
          gTr(i,j,k  ,iPIC ) = gTr(i,j,k  ,iPIC ) - flux/dzup
          flux = wC_sink*MAX(0.0, Ptr(i,j,k,bi,bj,iPOC))
          gTr(i,j,k  ,iPOC ) = gTr(i,j,k  ,iPOC ) - flux/dzup
          flux = wN_sink*MAX(0.0, Ptr(i,j,k,bi,bj,iPON))
          gTr(i,j,k  ,iPON ) = gTr(i,j,k  ,iPON ) - flux/dzup
          flux = wP_sink*MAX(0.0, Ptr(i,j,k,bi,bj,iPOP))
          gTr(i,j,k  ,iPOP ) = gTr(i,j,k  ,iPOP ) - flux/dzup
          flux = wSi_sink*MAX(0.0, Ptr(i,j,k,bi,bj,iPOSi))
          gTr(i,j,k  ,iPOSi) = gTr(i,j,k  ,iPOSi) - flux/dzup
          flux = wFe_sink*MAX(0.0, Ptr(i,j,k,bi,bj,iPOFe))
          gTr(i,j,k  ,iPOFe) = gTr(i,j,k  ,iPOFe) - flux/dzup
          DO l = 1, nplank
           flux = biosink(l)*MAX(0.0, Ptr(i,j,k,bi,bj,ic+l-1))
           gTr(i,j,k  ,ic+l-1 )=gTr(i,j,k  ,ic+l-1 ) - flux/dzup
#ifdef DARWIN_ALLOW_NQUOTA
           flux = biosink(l)*MAX(0.0, Ptr(i,j,k,bi,bj,in+l-1))
           gTr(i,j,k  ,in+l-1 )=gTr(i,j,k  ,in+l-1 ) - flux/dzup
#endif
#ifdef DARWIN_ALLOW_PQUOTA
           flux = biosink(l)*MAX(0.0, Ptr(i,j,k,bi,bj,ip+l-1))
           gTr(i,j,k  ,ip+l-1 )=gTr(i,j,k  ,ip+l-1 ) - flux/dzup
#endif
#ifdef DARWIN_ALLOW_SIQUOTA
           flux = biosink(l)*MAX(0.0, Ptr(i,j,k,bi,bj,isi+l-1))
           gTr(i,j,k  ,isi+l-1)=gTr(i,j,k  ,isi+l-1) - flux/dzup
#endif
#ifdef DARWIN_ALLOW_FEQUOTA
           flux = biosink(l)*MAX(0.0, Ptr(i,j,k,bi,bj,ife+l-1))
           gTr(i,j,k  ,ife+l-1)=gTr(i,j,k  ,ife+l-1) - flux/dzup
#endif
          ENDDO
#ifdef DARWIN_ALLOW_CHLQUOTA
          DO l = 1, nPhoto
           flux = biosink(l)*MAX(0.0, Ptr(i,j,k,bi,bj,iChl+l-1))
           gTr(i,j,k  ,iChl+l-1)=gTr(i,j,k  ,iChl+l-1)-flux/dzup
          ENDDO
#endif
         ENDIF
        ENDDO
       ENDDO
      ENDDO
#endif /* DARWIN_BOTTOM_SINK */

#endif /* ALLOW_DARWIN */

      RETURN
      END
