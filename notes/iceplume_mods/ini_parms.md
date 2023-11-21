## ini_parms

A small addition is in this script (around line 1540) after the mass2rUnit is calculated:
```
#ifndef ALLOW_ICEPLUME
C--   For backward compatibility, set temp_addMass and salt_addMass
C     to temp_EvPrRn and salt_EvPrRn if not set in parameter file "data"
      IF (temp_addMass .EQ. UNSET_RL) temp_addMass = temp_EvPrRn
      IF (salt_addMass .EQ. UNSET_RL) salt_addMass = salt_EvPrRn
#else
      temp_addMass = UNSET_RL
      salt_addMass = UNSET_RL
#endif
```
