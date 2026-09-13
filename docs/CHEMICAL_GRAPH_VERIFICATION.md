# Producer-owned covalent facts

The schema1 finite-chain inputs recorded atom connectivity but omitted bond
orders. A consumer cannot recover declared chemistry by guessing valence or
bond lengths. ChemicalBondGraph now owns a complete immutable ordered-atom
covalent unit; Structure's connectivity is derived from that one unit.
Pendant declares its internal orders (single by default, nitrile triple),
and the single scalar builder propagates them through backbone, attachment,
fragment and cap edges. Batch construction uses the same template owner.
Frame.to_structure consumes the explicitly requested polymer's declared graph
after element/connectivity validation; it does not infer orders from geometry.

Schema2 boundary manifests include these units, actual producer code hashes
and NumPy version. All nine XYZ byte payloads/canonical hashes and atom order
are unchanged. Existing Sarco imported schema1 manifests and completed native
journals remain untouched. This is a new chemistry receipt, not relabeling
the original calculations as having emitted these metadata.

137targeted chemistry/chain/fragment/forcefield/boundary tests pass. Sarco's
single finite_graph_identity decoder now accepts explicit atom-order bond
graphs as well as SMILES/AddHs order.31combined graph/stereo/host/preparation/
topology tests pass, and actual published CNEPO inputs still requalify using
their historical writer hashes. Native RDKit raw chirality tags refer to bond
insertion order; the explicit graph constructor canonicalizes that order.
Cross-encoding physical comparisons use CIP and ligand-ordered coordinate
parity, not raw tag equality. No duplicate stereochemistry algorithm exists.

Actual old-input/new-producer metadata bridging and whole-geometry checks
pass for all nine sampled-return-qualified Sarco references. PVDF/VDCN have
no stereocenters. AN source-order atoms15/21/27 at5/7/9 lengths are R and
remain so at the selected endpoints. No minimum/Hessian, field response,
crystal packing or physical material gate follows from preserved chemistry.

## Full-suite baseline diagnosis

The running full suite exposed three frozen-energy literal comparisons:
PE, alpha and gamma in test_mechanics.py. An isolated detached baseline at
b600c33 reproduces all three failures before this change. Baseline/current
zero-field and fixed-field energies for PE/alpha/beta/gamma are exactly equal
as binary float hex values. On this runtime, literal mismatches are
PE -4.143020937018385 vs recorded -4.143020937018386;
alpha27.71522759686365 vs27.715227596863407;
gamma27.706195941151186 vs27.70619594115115.
No test literal, numerical tolerance, physical gate or kernel was changed to
hide these preexisting failures. Full-suite final counts are recorded after
its own process completes; the full suite is not claimed green.

The completed first full-suite run reports580passed,5skipped and exactly
the three baseline failures in331.13s. A subsequent immutable-constructor
regression snapshots caller-owned pendant atom/bond/order collections;
its fresh verification is recorded separately, not retroactively assigned
to that earlier full-suite run.

The fresh full suite after that constructor snapshot change completes with
581passed,5skipped and the same three baseline literal failures in329.78s.
No new failure was introduced. Only the compact baseline diagnosis remains;
the temporary detached source checkout was removed, recoverable at b600c33.
