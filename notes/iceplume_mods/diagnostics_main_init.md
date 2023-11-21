# diagnostics_main_init.F

Add a new block for addMass around line 178:
```
#ifdef ALLOW_ADDFLUID
      diagName  = 'addMass '
      diagTitle = 'Source (<0: sink) of fluid in the domain interior'
      diagUnits = 'kg/s            '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
c
      diagName  = 'gTaddMas'
      diagTitle = 'Temperature tendency of addMass, >0 inc ocean theta'
      diagUnits = 'degC/s          '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
c
      diagName  = 'gSaddMas'
      diagTitle = 'Salinity tendency of addMass, >0 inc ocean salt'
      diagUnits = 'g/kg/s          '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
      diagName  = 'tAddMass'
c      diagTitle = 'Temperature of AddMass in apply_forcing'
      diagTitle = 'Temperature times AddMass in apply_forcing'
      diagUnits = 'C.kg/s.m3/kg/m3 '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
c
      diagName  = 'sAddMass'
c      diagTitle = 'Salinity of AddMass in apply_forcing   '
      diagTitle = 'Salinity times AddMass in apply_forcing   '
      diagUnits = 'g/kg.kg/s.m3/kg/m3'
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I     diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* ALLOW_ADDFLUID */
```
