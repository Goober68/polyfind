# Clamped-ion piezoelectric columns of beta-PVDF (Polyfind model)

A reproducible record of the consumer side of the clamped-ion comparison in
`docs/REFERENCE_DATA_REQUEST.md` (2026-09-13), produced by
`examples/clamped_ion_columns.py`. It is a model quantity on the fitted
potential (`pvdf-dft-valence-flux-born`, induced dipoles, Ewald), not a fit,
not a material coefficient and not a same-Hamiltonian quantity. The producer's
Berry-phase values are not reproduced here and carry their own labels.

`beta_pvdf_clamped_ion_columns.json` holds the exact reference cell and packer
parameters, the placed geometry in the producer's frame (also as
`beta_pvdf_reference_producer_frame.xyz`, SHA-256 `c056e4901828bc107967a1ca2e5476457fbd139bdb6d60b84d7297cc8d2725e4`),
the frame map, the charge-only / induced / total polarization, both definitions
of every expressible column (dipole per reference volume, and Vanderbilt's
proper tensor with the `- delta_ij P_k` term) at two amplitudes, the shear
convention, and the SHA-256 of every source file that defines the calculation.
The `yz` and `xz` shears involve the chain axis and are not expressible in this
cell parametrisation; they are recorded as such.

Regenerate with:

```text
python examples/clamped_ion_columns.py
```
