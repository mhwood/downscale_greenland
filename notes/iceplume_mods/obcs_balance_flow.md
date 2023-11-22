# obcs_balance_flow.md

This modification is the only change made by me to ensure that mass is note accumulated in the domain resulting from the plume (or any addMass-derived mass flux for that matter - I'm a bit surprised no one else ran into this issue).

First, a new variable needs to be defined under the shelfIceNetMassFlux variable (around line 58):
```
      _RL addMassMassFlux
```
Then initialize the addMassMassFlux before the OBCSbalance loop (around line 290):
```
       addMassMassFlux = 0. _d 0
```

Then, a loop is added so that the addMass flux is accounted for in the sum of the imbalance (around line 320):
```
#ifdef ALLOW_ADDFLUID
      IF ( selectAddFluid.NE.0 ) THEN
        DO bj=myByLo(myThid),myByHi(myThid)
        DO bi=myBxLo(myThid),myBxHi(myThid)
           tileFlow(bi,bj) = 0.
           DO j=1,sNy
           DO i=1,sNx
           DO k=1,Nr 
             tileFlow(bi,bj) = tileFlow(bi,bj)
     &          + addMass(i,j,k,bi,bj) * maskInC(i,j,bi,bj)
           ENDDO
           ENDDO
           ENDDO
        ENDDO
        ENDDO
        CALL GLOBAL_SUM_TILE_RL( tileFlow, addMassMassFlux, myThid )
        IF ( debugLevel.GE.debLevC ) THEN
          WRITE(msgBuf,'(A,I9,A,1P1E16.8)') 'OBCS_balance (it=',
     &       myIter, ' ) correct for addMassMassFlux:',
     &       addMassMassFlux
          CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &       SQUEEZE_RIGHT, myThid )
         ENDIF
      ENDIF
#endif /* ALLOW_ADDFLUID */
```

Then, the component from addMass is added after the shelfice component (around line 377):
```
#ifdef ALLOW_ADDFLUID
         IF ( selectAddFluid.NE.0 )
     &        inFlow = inFlow + addMassMassFlux*mass2rUnit
#endif
```
