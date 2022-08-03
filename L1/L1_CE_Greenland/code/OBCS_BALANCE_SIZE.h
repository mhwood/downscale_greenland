C------------------------------------------------------------------------------|
C                           OBCS_BALANCE_SIZE.h
C------------------------------------------------------------------------------|
C The determines the number of timesteps that can be recorded
C This determines the maximum size of the vector which records flux imbalances
C The first N values are averaged when OB_balancePeriod > 0 
C The length needs to be a minimum of ceil(OB_balancePeriod/timestep)

      INTEGER, PARAMETER :: nOBCS_fluxImbalance_rec = 6720

