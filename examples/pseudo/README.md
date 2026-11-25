# Three pseudopotential files `examples/chain_triage.py` uses

The other files of its set-up (H, C, F, N, Cl) are fetched from the SSSP 1.3.0 PBE precision
library by `examples/gpu_jobs/dft_cells.sh` and `examples/mi100_setup.sh`, with their checksums
checked. These three are not in that library, so they are kept here: a number computed on another
machine is then computed with the same files.

| File | What | From |
|---|---|---|
| `O.pbe-kjpaw.UPF` | oxygen, PAW, PBE: relaxations and energies | Quantum ESPRESSO's own pseudopotential set (Debian package `quantum-espresso-data` 6.7), GPL as that distribution |
| `C.pbe-mt_gipaw.UPF` | carbon, norm-conserving, PBE: Berry-phase runs only | the same set |
| `O.pbe-tm-triage.UPF` | oxygen, norm-conserving, PBE: Berry-phase runs only | made here by `ld1.x` 7.5 from `examples/chain_triage_oxygen.in` |

Why two sets: the Berry phase of `pw.x` 7.5 is right only with norm-conserving files (a third of a
quantum off on a lone PVDF chain with the PAW carbon), so polarisation runs take the
norm-conserving carbon and oxygen on the geometry the PAW files relaxed.

sha256:

    2b6ba2a535fad72e0dd8046da3f14876a13829f33e398447fc366012bd7b8999  O.pbe-kjpaw.UPF
    a30eca0fbeec806e3015e87dd58d945146f7e10b3337380053d9bf4089e71504  C.pbe-mt_gipaw.UPF
    c8d54f0c9425650b32bb4f1e9e56b009839ab50f524176a5122d8a952ddac3ba  O.pbe-tm-triage.UPF
