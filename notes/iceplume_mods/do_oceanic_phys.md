# do_oceanic_phys.md

There are only two changes to this script.

First, a note is added for ICEPLUME_CALC in the pipeline notes just after ICEFRONT_THERMODYNAMICS:
```
C       |
C       |-- ICEPLUME_CALC
```

Second, the ICEPLUME_CALC function is added just after the ICEFRONT function (around line 477):
```
#ifdef ALLOW_ICEPLUME
      IF ( useICEPLUME .AND. fluidIsWater ) THEN
C     Calculate plume due to subglacial runoff, as well as 
C     temperature and virtual salt flux at ice-ocean interface
       CALL TIMER_START('ICEPLUME_CALC [DO_OCEANIC_PHYS]',
     &       myThid)
       CALL ICEPLUME_CALC( myTime, myIter, myThid )
       CALL TIMER_STOP( 'ICEPLUME_CALC [DO_OCEANIC_PHYS]',
     &      myThid)
      ENDIF
#endif /* ALLOW_ICEPLUME */
```
