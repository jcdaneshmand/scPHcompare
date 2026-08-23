# MV8-SA baseline-helper recovery prefreeze

The exact MV8-S execution at `218af286...` stopped on its first child before
cell selection because the child omitted the helper defining preserved RNG state.
It ran 32.7 seconds, stayed below both caps, published no baseline, and ran no PH.

This amendment adds only the missing helper import and authorizes one complete
replacement execution in fresh roots after commit. The failed roots are immutable.
All scientific axes, transforms, PH definitions, fallback/cap rules, landscape
definition, and downstream firewalls are unchanged.
