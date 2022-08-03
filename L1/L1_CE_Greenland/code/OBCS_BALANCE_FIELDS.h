#ifdef ALLOW_OBCS

CBOP
C     !ROUTINE: OBCS_BALANCE_FIELDS.h
C     !INTERFACE:
C     #include "OBCS_BALANCE_FIELDS.h"

C     !DESCRIPTION:
C     *==========================================================*
C     | OBCS_BALANCE_FIELDS.h
C     | o Header file containing flux imbalances
C     *==========================================================*
CEOP

#ifdef ALLOW_OBCS_SMOOTH_BALANCE
      COMMON /OBCS_BALANCE_FIELDS/
     & OBCS_fluxImbalances
      _RL OBCS_fluxImbalances(nOBCS_fluxImbalance_rec)
#endif /* ALLOW_OBCS_SMOOTH_BALANCE */


#endif /* ALLOW_OBCS */
