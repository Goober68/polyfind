# Finite-chain boundary-sensitivity starts

These are chemistry-owned inputs for Sarco stage 2, not relaxed structures or
predicted field responses. The 5/7/9-monomer ladder adds VDF host units
symmetrically around a fixed center. AN and VDCN therefore contain exactly one
central defect at every length; PVDF contains none.

Each chain is neutral, H-capped, connected, and initially all-trans. Coordinates
use `z` for the backbone end-to-end direction, `y` toward the central substituted
unit's non-H pendants, and `x` to complete the right-handed transverse frame.
`manifest.json` owns the sequence, stereochemistry statement, complete bond graph,
atom and backbone indices, axes, canonical-LF SHA-256, and component validation
for every XYZ. Line-ending changes between Windows and Linux do not change identity.

The intended Sarco protocol independently varies:

- length: 5, 7, then 9 monomers;
- the number/registry of neighboring parallel chains;
- endpoint-clamped versus rotation-permitted boundaries;
- field along `x`, `y`, and `z`.

The size gate requires the same qualitative mechanism and no more than 10% change
in the extrapolated observable across both 5->7 and 7->9 increments. These
coordinates do not themselves clear that gate.

CNEPO is deliberately absent. Its epoxide bridges two backbone atoms and changes
the backbone bond graph; Polyfind's pendant-fragment model cannot represent it.
Treating the oxygen as an ordinary pendant would create the wrong molecule.

Regenerate with:

```text
PYTHONPATH=src python examples/boundary_sensitivity_starts.py
```
