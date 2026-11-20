# Resume: polyfind, state as of 2026-09-13 (afternoon)

**STANDING OBJECTIVE (PI, 2026-11-14; also in `CLAUDE.md`).** "strain is NOT the sole metric by which to gauge. What we want
is the to be able to evaluate candidates for the optimal point on the strain/tanD/freq curve. For example: a high frequency,
low loss moderate strain polymer will outperform our CFE." Every candidate is reported and ranked on strain, tan delta and
frequency together; one with less strain than the CFE terpolymer is not out. The window reading of section 3 of
`docs/CANDIDATE_MATERIALS.md` is strain only and is a partial result by this rule; the columns that carry all three are the
frequency maps' (`tan delta per % strain`, `efficiency`).

**Final status of the 14-hour window, 2026-11-07 19:06 UTC.** Closed today: the stage-7 problem solved by bounded
least squares (`pvdf-dft-valence-belt3`, 5.18-5.19); the sensitivity, screen and basin analyses (5.15, 5.17, 5.20); the
Gamma-route relaxed-ion tensor (5.16); the pendant-group finding and its fix (5.22: the rigid refinement cannot tilt the
CH2/CF2 groups as DFT's alpha does, the gap row carried 0.23 kcal/mol per monomer of that strain, `polyfind.pendants`
relaxes the groups as rigid tilts, the refit is `pvdf-dft-valence-belt4`, the sixth prediction set
`e1/e2_polyfind_belt4.json` carries the correction, the copolymer model reorders to TGTG' lowest); the soft C-C stretch
constant quantified on the Gamma route (beta unmoved, alpha's polar eps_0 3.47 -> 2.93 floored); stage 8 (shape term,
5.21) and stage 9 (tilts as refinement variables, 5.23) designed. Open: the lateral coupling (C_12 0.87, the one residual
past a tolerance at belt4), waiting on Dipole's rotated-setting cell for its anchor; Dipole's alpha and Gamma successors
(launched 05:10 UTC, no result by 19:06) -- when they land: `examples/dft_cell_geometry.py` on the alpha record for the
backbone angles and the pendant tilts (prediction CF2 wag 2.0-2.5, rock 1.1-1.4, twist 0.6-0.7 deg), the Gamma DFPT
against 2.23/4.61/2.73 and the Born charges, then re-evaluate belt4 and refit only if a residual leaves 1.5 tolerances
(task #16). Suite 724 passed, 5 skipped. Branch `claude/polymeric-stable-arrangements-uh06b1`, head 7f68362 before this status commit.

**POLYFIND'S CELL MAY NOW SHEAR ALONG ITS CHAINS, AND IT DOES NOT FOLLOW DFT: POLYFIND CHARGES THREE TIMES DFT FOR
MOVING CHAINS ALONG EACH OTHER, 2026-11-18 (PI: "Yes" to "give polyfind's cell the two missing shear freedoms and
re-polish the candidates").** Text: `docs/MLIP_SECOND_OPINION.md` 3.9 (it corrects 3.8's "polyfind's energies are
right": true at one registry, not along it).
(1) THE FREEDOM, DONE AND EXACT. `polyfind.mixedcell.MixedModel(shear=True)`: two more variables at the END of the
    vector (`index("shear")`), what a and b carry along the chains, in A; `lattice_rows(A, B, gamma, c, sa, sb)`;
    `MixedState.shear` (None = upright; `from_vector(..., shear=True)`); `tile_state` scales it; `reduce_shear` takes
    whole repeats off. Upright vectors and every stored number are untouched (bit for bit at no shear). Known answer:
    a cell sheared by half a repeat = the upright cell twice as wide, 1e-14. `tests/test_mixedcell.py` (30 shear tests).
    TWO TRAPS: (a) many stored cells have a symmetry that makes no shear a stationary point, so a relaxation must
    start off it (`shear_polish.py --kick`, 0.05 A); (b) SLSQP's own stopping test goes on failing at the minimum of a
    14-monomer cell (1300 evaluations seen): use `relax(..., patience=60, patience_tol=1e-6)`.
(2) THE CANDIDATE (`vdf5-rvdf1-6` all-trans; `examples/shear_registry.py candidate`,
    `docs/predictions/shear_registry_belt4.json`): DFT -0.35 kcal/mol per monomer at shear (0.30, -0.90) monomers;
    polyfind with the shear free -0.03 at (0.026, -0.078), the same direction, a tenth of the way. DFT's relaxed cell
    fits the mixed cell's variables to 0.015 A; polyfind there is +0.43 (shear held, the rest relaxed), +1.33 as fitted.
(3) WHY, MEASURED (`shear_registry.py frames` + `mlip_frames.py label` + `collect`; 35 labelled frames in
    `docs/predictions/mlip/shear_registry_frames_pbe_d3bj.extxyz`; `tests/test_shear_registry.py`): plain PVDF,
    twelve-atom cells sheared rigidly, two packings, polyfind against PBE-D3(BJ) single points on the same atoms:
    polyfind is 2.4 to 3.5 times DFT (median 2.9 over 32 cells; half a monomer along a in beta: 2.34 against 0.69).
    THE REGISTERED E5 TEST, run here in the same DFT (Dipole never delivered it; the forces on Dipole's beta cell are
    zero in ours): polyfind 0.172 / 0.656 / 1.251 / 1.533 (one sublattice slid c/8 .. c/2) against DFT 0.066 / 0.230 /
    0.396 / 0.467, and one chain of four 0.671 / 1.585 against 0.232 / 0.469. Sign, order and anisotropy hold; the
    factor-1.5 acceptance is missed (2.6 to 3.4). It is all in polyfind's repulsion and dispersion pair terms (its
    electrostatics lower the sheared cell). Handles (`shear_registry.py handles`): the H diameter, the H-F pair's
    diameter, the belt wall's depth, the F diameter, which are the lattice's handles too; as a linear estimate no
    moderate change of them lowers the cost and holds the lattice (a softer wall for H's contacts is the likely need).
(4) WHAT STANDS: energies at one registry (beta-alpha gap, the candidate's margin with both cells upright, C_11 C_22
    C_12). WHAT DOES NOT, UNTIL REDONE: every barrier of a chain or a kink moving along its neighbours (the 3.15
    chain-slide barrier, the kink steps behind every relaxation frequency and tan delta: overstated, by how much not
    computed, so polyfind's rates are too slow), the shear constants along the chains, and the margins of cells with
    pins (a pin's registry with the next chains' pins is worth about the window's width; DFT moves it, polyfind cannot).
(5) POLYFIND'S OWN RE-POLISH WITH THE SHEAR FREE (`examples/shear_polish.py`, resumable, `--minutes`;
    `docs/predictions/sheared_margins_belt4.json`): DONE, all 35 rows, 168 cells,
    every one stationary to 1e-4, rigid (tilts NOT redone). The all-trans cells stay upright (mean gain 0.004, largest
    0.029 = the reversed-unit candidate's, which DFT lowers by 0.35); the helix cells move (mean 0.038, largest 0.294,
    38 of 133 over 0.05, half of those near half a monomer of shear). Rigid margins: -0.028 to +0.265, mean +0.044;
    9 rows move by more than 0.06, 4 by more than the window's 0.12 (`vdf9-trfe3-an1-rvdf1-14` +0.265, `vdf6-trfe4`
    +0.183, `vdf6-trfe3-cfe1` +0.141, `vdf8-trfe4-an1-rvdf1-14` +0.124); the CFE terpolymer +0.111; the lowest helix
    conformer changes in 16 rows. NOT A CORRECTION TOWARD DFT (polyfind -0.03 where DFT is -0.33, and on the other
    side): it shows that the registry a cell is allowed moves a margin by the window's width. The polished margin table
    (`polished_margins_belt4.json`, with tilts) is NOT replaced and the ranking by margin is not settled by either.
(6) DFT WORKING DIRECTORIES (outside the repository, survive restarts): `/opt/dft/shear2` (the 28 sheared cells and
    the seven E5 cells, start and labelled), `/opt/dft/shear`, `/opt/dft/e5` (the first labelling; the same frames).
(7) THE REFIT, IN PROGRESS (PI 2026-11-18: "Yes" to "Refit polyfind's non-bonded wall with registry in the target
    set, probably with a softer wall for hydrogen's contacts?"). After it: the slide and kink barriers and the
    frequency maps again; margins of cells with pins as a spread over registries, not one upright cell.
    (a) DONE: stage 11 of the crystal fit (8 and 9 are CRYSTAL_STAGE.md's unbuilt designs) (`crystalfit.REGISTRY`, `registry_costs`, `evaluate(registry=True)`,
        `examples/crystal_stage.py --stage 11 [--wide]`): eight registry rows beside the sixteen anchors, tolerance a
        fifth of DFT's value. belt4 scores 699 on it (11 on the old objective). TRAP: inside `applied()` build packers
        as `pack_mod.CrystalPacker(...)`, never from a name imported beforehand (the block rebinds the class on the
        module; an early import gives the polymer's own charges and the default torsion term, silently).
    (b) DONE, AND IT IS THE NEGATIVE ANSWER: with the present form the fit brings the registry cost to 1.0-1.3 of DFT's
        and the lateral cell contracts 2-3% to pay (objective 43 inside the old bounds, 46 with the diameters'
        bounds opened; `docs/predictions/crystal_stage_v11_registry[_wide].json`, both stopped by hand at 186
        evaluations, still falling slowly, NOT presets). The form cannot hold both.
    (c) DONE: the data a wall is fitted to on one footing. `examples/interchain_frames.py frames --out
        /opt/dft/interchain [--packings ...]`: plain PVDF, chains RIGID, lateral edges scaled -6..+6% singly and
        together, shears; packings `beta` (27 cells), `candidate` (27), `alpha` = TGTG' `vdf2` (27), and `hh` = the
        chain with every second unit reversed, -CH2-CF2-CF2-CH2- (`vdf1-rvdf1-2`, 37 cells: the only ones with
        fluorine against fluorine and hydrogen against hydrogen between chains). Label: `python
        examples/mlip_frames.py label --frames /opt/dft/interchain --pw /opt/mm/qe/bin/pw.x --pseudo /opt/pseudo
        --np 4 --minutes 100` (resumable; 50 s per 12-atom cell, 2.3 min per 24-atom one). ALL LABELLED (beta 27, candidate 27,
        alpha 27, hh 37, cand6 4) and collected (`interchain_frames.py collect --frames /opt/dft/interchain`).
        WHAT BETA SHOWS: polyfind (belt4) is right on the stretched side and 1.2-1.4x too stiff at -6%; its registry
        cost over DFT's is 2.2x at a 4% wider cell, 3x at the cell, 5x at a 4% narrower one: too steep close in.
    (d) DONE: `examples/interchain_wall.py pairs --frames /opt/dft/interchain --registry /opt/dft/shear2` (the
        between-chain distances of every cell by element pair; the present form on them is polyfind's own pair
        energy to 4e-14; RERUN IT after new labels, it reads the DFT energies) and `... fit --form
        12-6|n-6|exp-6|b14-7 --free ...`. On 52 cells, depths held: 12-6 0.12 kcal/mol per monomer rms, n-6 0.04,
        exp-6 0.014, buffered 14-7 with its buffer free 0.007 and nothing at a bound. Variables `x_FF, D_FF, x_HH,
        D_HH` give the like pairs their own terms (for the hh cells).
    (e) DONE: THE NEW WALL IN THE KERNELS, AS AN OPTION. `pack.wall_b147`, `CrystalPacker(wall=("b14-7", delta))`
        (two tables = depth and inverse diameter; NEGATIVE depth = a wall alone = the belt sites, which stay),
        mixedcell's batched kernel, supercell's sum, `FFParameters.wall` forwarded by `applied()`. `wall=None` is
        12-6, untouched. Checked against the standalone trial to 4e-11; `tests/test_wall.py` (11). Its attraction
        falls as r^-7 (an effective form inside the cutoff). The tabulated screen stays 12-6 (it only picks starts).
    (f) DONE, INTERMEDIATE: stage 12 = stage 11 with the wall, its buffer the 14th variable
        (`crystal_stage.py --stage 12 --geometry pbe --method trf --pendant-tilts --x0 <belt4 file> --set
        x_C=...`; `docs/predictions/crystal_stage_v12_wall.json`, stopped by hand at 240 evaluations, NOT a
        preset): objective 20.0 (11.1 of it ridge), gap -0.730, beta cell within 0.3%, C_11 C_22 C_12 22.0 22.6 6.1,
        all 8 registry rows within 7%, alpha b 1.7% short (the one row past a tolerance). HELD OUT: the 18 other
        sheared plain cells 0.97-1.08 of DFT (belt4 2.4-3.5).
    (g) THE REVERSED CHAIN (`hh`, 37 cells, ALL LABELLED): a buffered wall fitted with the plain cells leaves 0.154
        kcal/mol per monomer rms on them (plain groups 0.01-0.04); F-F and H-H pair terms of their own DO NOT help
        (0.153, three bounds). The misses are long shears where polyfind's electrostatics alone exceed DFT's whole
        energy. `interchain_wall.py fit` can now move the electrostatics: `es` (a multiplier, diagnostic: hh alone
        asks 0.38, plain 0.92), `s_CH`, `s_CF` (the two groups' charges), `d_CC` (charge along each CH2-CF2 backbone
        bond), every cell's electrostatic energy stored as its exact quadratic form in them (`K`, induction in it),
        and the extra sources of `extra_sources` to first order (`K1`, a plain Ewald pair matrix `ewald_green`).
        ENERGIES CANNOT SET CHARGES: such fits end on bounds (CF2 charges to zero).
    (h) THE POTENTIAL, MEASURED (`examples/chain_potential.py frames|dft|collect|fit`, working directory
        `/opt/dft/chainpot`, DFT 3 to 6 minutes a chain on one core; stored samples
        `docs/predictions/chain_potential_pbe_samples.npz`, fit `chain_potential_belt4.json`,
        `tests/test_chain_potential.py`): one chain alone in a 16 A box, DFT's electrostatic potential (pp.x
        plot_num 11) 2.2-5 A from the nuclei against polyfind's charges + induced dipoles. Error 0.41 / 0.38 / 0.42
        kcal/mol per e rms around plain all-trans / plain TGTG' / the reversed chain = 3.6% / 5.5% / 47% of the
        potential (the reversed chain has no dipole). No atom charges repair it, no single term on F; A TERM ON
        HYDROGEN DOES on all three (a dipole or quadrupole along C-H, or a charge site on the bond): 0.18 / 0.13 /
        0.23. Text: `docs/CRYSTAL_STAGE.md` 5.26.
    (i) TWO OPTIONAL CHARGE TERMS IN THE KERNELS, OFF BY DEFAULT, NO PRESET USES THEM (`tests/test_bond_sites.py`):
        `polymers.with_bond_sites(polymer)` puts a site `M` on every C-H bond (0.65 of the way out; charge only: no
        LJ, no polarizability, no dipole, excluded and scaled as its hydrogen, not counted in the chain's frame, so
        a sited chain sits exactly where the plain one does and old stored frames still place it);
        `FFParameters.bond_sites=(("H", q),)` charges it (and `applied()` REFUSES a chain without sites);
        `FFParameters.polymer(p)` builds what a set needs; examples: `candidate_screen.SITES_FOR = params` makes
        `make_polymer` build sited chains. `FFParameters.backbone_pull=(("F", e),)`: each F draws e onto its carbon
        from each backbone neighbour (= the trial's `d_CC` at -2e; tells CF2-beside-CH2 from CF2-beside-CF2).
        With charges from the potential alone and only the wall refitted: hh cells 0.154 -> 0.128 (site) -> 0.089
        (site + pull 0.0155: dCH -0.049, dCF +0.118, site -0.318). TRAP: `applied()` blocks do not nest (the outer
        set's charges win).
    (j) THE CANDIDATE'S OWN CELL, RIGID: THE GAIN IS ALL BETWEEN CHAINS (`interchain_frames.py` packing `cand6`:
        `vdf5-rvdf1-6` all-trans polished cell, own chains rigid, sheared along DFT's line; 72 atoms, 21 minutes of
        DFT a cell; ALL FOUR LABELLED). DFT ON THE RIGID CELLS (quarter, half, end of the line to (0.30, -0.90)
        monomers): -0.359, -0.412, -0.338 kcal/mol per monomer (relaxed, at the end: -0.346). REGISTERED BEFORE
        THEM (`docs/predictions/interchain_cand6_predictions.json`): belt4 +0.10, +0.70, +1.58; buffered wall,
        present charges -0.21, -0.14, +0.15; with the bond site -0.24, -0.21, +0.03; with site and pull -0.25,
        -0.23, -0.04. SO: the wall takes 1.4 of the 1.9 at the end, the charges 0.2; 0.10 / 0.18 / 0.30 is left, in
        F..F at 2.52 A (8 a cell, level reversed units) traded for H..F at 2.42.
        ALL 150 RIGID CELLS ARE LABELLED AND STORED (`docs/predictions/interchain_belt4.json`, frames beside it,
        `tests/test_interchain.py`), AND THE TRIALS ON THEM (`docs/predictions/interchain_wall_trials.json`, 11
        trials; hh rms / plain groups / the candidate's line): buffered wall, present charges 0.17 / 0.03-0.09 /
        -0.28 -0.26 -0.01; site + pull at the potential's strengths 0.09 / 0.01-0.06 / -0.29 -0.30 -0.13; with the
        F-F pair's own diameter and depth (`--free ...,x_FF,D_FF`) 0.10 / 0.01-0.03 / -0.33 -0.38 -0.20; that with
        `--weight cand6=4` 0.13 / 0.02-0.05 / -0.36 -0.44 -0.29; the same weighting WITHOUT the charge terms hh
        0.25. Rebuild the trial's tables after any new label: `interchain_wall.py pairs --frames
        /opt/dft/interchain --registry /opt/dft/shear2 [--charges M=-0.3 --out /opt/dft/interchain/pairs_sites.npy]`.
    (k) THE CONFORMER SIDE, STARTED: `ValenceDesign.lattice = <crystal-stage FFParameters>` makes the 467 conformers
        (`POLYFIND_TRAINSET=/home/user/sarco/materials/gpu_bundle/trainset_v1_566.xyz`) see the lattice's
        non-bonded terms; `fit_valence(designs=(design, held), free=range(12), x0=F.VAL_FITTED_X)` refits the
        torsions under them in 9 s. Energy rms train/held-out: valence fit 1.02/1.36; belt4 1.31/1.72, torsions
        refitted 1.21/1.66; stage 12 1.68/1.97, refitted 1.49/1.84. NOT ACTED ON: no preset carries refitted
        torsions. TWO THINGS THIS EXPOSED, BOTH OLDER THAN THE REFIT: the crystal stage never refitted the torsions
        to its own non-bonded terms, and the lattice side uses ONE torsion triple (the CF2-CH2 type's) for every
        backbone bond of every copolymer (`reference_ff_parameters`), though the valence fit has four types. The
        conformer side does not know the bond sites or the pull yet (SimpleFF leaves both out).
    (l) THE CONFORMER SIDE KNOWS THE TWO CHARGE TERMS (`ValenceDesign.lattice` with `bond_sites`/`backbone_pull`;
        off-nucleus F charge then projected as the packer does, at `fitting.PROJECT_BOND`). Torsions refitted on
        the 467 conformers under each set (train / held-out energy rms): valence fit 1.02 / 1.36; belt4 1.21 /
        1.66; stage-12 wall 1.49 / 1.84; that wall with the potential's charges 1.30 / 1.63 and the best force rms
        of all (4.600 / 5.207). BUT the conformer-refitted torsions leave the alpha-beta gap at -0.24 (DFT -0.726).
        Holding the one-fold CF2-CH2 term at -1.6 / -2.0 / -2.4 and refitting the rest: conformers 1.303 / 1.312 /
        1.328 (held-out improves), gap -0.33 / -0.43 / -0.55, alpha b 9.53 / 9.59 / 9.66 (DFT 9.636). THE
        CONFORMERS LEAVE THAT DIRECTION OPEN AND THE CRYSTAL SETS IT. Text: `docs/CRYSTAL_STAGE.md` 5.27.
    (m) STAGE 13 IS ONE FIT: 21 variables (`crystalfit.VARIABLES13` = stage 12's 14 + `q_M`, `pull_F` + `V1 V2 V3`
        (the lattice's torsion triple, no ridge) + `x_FF`, `D_FF`), rows = stage 12's 24 + `examples/stage13_rows.py`
        (`Rows()(params)`: 117 rigid cells in the kernels, 3 potentials, 2 conformer rows; needs `POLYFIND_TRAINSET`).
        Command: `POLYFIND_TRAINSET=/home/user/sarco/materials/gpu_bundle/trainset_v1_566.xyz python -W ignore
        examples/crystal_stage.py --stage 13 --geometry pbe --method trf --pendant-tilts --rows --x0 <result or
        start json> --maxfev 10 --out /opt/dft/crystal_stage_v13.json` (16 s an evaluation, 22 a Jacobian; the
        --out file is the checkpoint; `python examples/stage13_rows.py <fit json>` prints the outside rows).
        FIRST RUN DONE (`docs/predictions/crystal_stage_v13_run1.json`, objective 520 -> 81, NOT converged, NOT a
        preset): lattice within 0.7%, gap -0.772, C 21.6 22.0 7.2, rigid cells 0.05-0.06 plain, 0.12 reversed
        chain, candidate's line -0.30 -0.35 -0.24, potentials 0.13-0.22, conformers 1.32; plain registry rows now
        0.65-0.78 of DFT (traded). SECOND RUN DONE (`docs/predictions/crystal_stage_v13_run2.json`, 81 -> 62, ten
        more iterations, still falling 1-3 an iteration, NOT a preset): gap -0.730, C 23.1 22.8 7.1, registry rows
        0.76-0.86 of DFT, rigid cells 0.03-0.04 plain and 0.13 reversed chain (unmoved: the terms do not reach it),
        candidate's line -0.32 -0.38 -0.26, conformers 1.29; paid with beta a 1.0% short. Not tried in the kernels: the trial's other basin (hydrogen smaller, pairs on
        their own terms).
    (n) A CHECK NO FIT SAW (`examples/cohesion.py`, `docs/predictions/cohesion_rigid_cells.json`): cell minus its own
        chains alone, DFT beta -7.82 / alpha -5.07 / hh -3.11 kcal/mol per monomer; belt4 -5.75 / -3.12 / -1.62; the
        wall -6.85 / -4.14 / -2.46; stage 13 run 1 -6.83 / -4.19 / -2.20; alpha minus beta right to 0.1 in all.
        METHOD NOTE: a DFT single point on polyfind's rigid geometry carries strain that differs between CONFORMERS
        by 0.3-0.5 kcal/mol per monomer (rigid-geometry alpha minus beta is -0.16 in DFT, -0.726 relaxed): compare
        only within a packing, or a cell with its own chains.
    (o) THE ACCEPTANCE ON `vdf5-rvdf1-6` HAS TWO OUTCOMES (`examples/candidate_acceptance.py --fit <json>
        [--torsions docs/predictions/valence_torsions_stage13.json]`, 7 minutes; stored
        `docs/predictions/candidate_acceptance_stage13_run1[_typed].json`): PASSES what the refit was for (all-trans
        shear gain -0.309, DFT -0.346, belt4 -0.03; helix none) and LOSES what belt4 had: all-trans minus DFT's
        helix (TG'TG) -0.135 upright / -0.444 relaxed (DFT +0.070 / -0.261; belt4 +0.064 upright), and GTG'T minus
        TG'TG -0.158 (DFT +0.197, belt4 +0.12); typed torsions from the conformer refit make that -0.46. The shear's
        part of the margin is right (-0.31 against DFT's -0.33); the LEVEL inside a chain with a reversed unit is
        off by 0.2-0.36. Split at the stored cells: the within-chain non-bonded energy of GTG'T against TG'TG moved
        by -0.19 between belt4 and run 1 (between chains +0.05). Nothing in the fit holds it; the conformer rows are
        too coarse (1.3 rms on an oligomer = 0.3 per monomer).
    (p) THE INSIDE OF A CHAIN HAS ITS ANCHORS NOW, AND THEY CHANGE THE READING (`docs/CRYSTAL_STAGE.md` 5.29).
        * A CELL SPLIT IN TWO: `examples/chain_split.py frames|dft|polyfind` takes a chain out of a DFT-relaxed
          cell alone into a 16 A box (one more single point; `mlip_frames.py label --pools 1` for a 36-atom chain,
          else out of memory); `docs/predictions/chain_split_pbe_d3bj.json`. Plain alpha-beta: cell -0.726 =
          inside -3.575 + between +2.850 (belt4 -3.62/+2.88; stage 13 -3.90/+3.16). Candidate, all-trans minus
          TG'TG: +0.070 = +1.060 - 0.990 (belt4 +1.71 - 1.64: two errors cancelling); GTG'T minus TG'TG: +0.197 =
          +0.415 - 0.218 (stage 13: -0.19 + 0.04).
        * A LONE CHAIN RELAXED FOR ONE MORE DFT POINT: `examples/mlip_train.py relax --dft-gradient` (our potential
          tilted to DFT's forces and stress at a labelled frame = a Newton step; right to 0.01 kcal/mol on 6- and
          12-atom chains; on a 36-atom chain it RUNS AWAY along soft turns: use `--max-move 0.3` and label again).
        * POLYFIND'S LONE HELIX FALLS SHUT: gauche angle alone DFT 56 (60 in the cell), belt4 41, stage 13 21-25.
          Stage 13's torsion triple (-1.9, 1.4, 0.9) was paying for an error between chains; fitted to the lone
          chain's line alone it is (0.77, -0.44, 0.66) with the angle stiffness 167 (not 264), and then the gap is
          -0.19 unless the rest moves.
        * THE LINE: `examples/chain_stretch.py scan|collect|polyfind`, `docs/predictions/chain_stretch_pbe_d3bj.json`
          (lone helix held at 4.35 to 5.18 A, lone all-trans chain at 2.525 to 2.66 A; tension along the chain in
          the file; chain stiffness about 130 GPa helix, 290 GPa all-trans). ELEVEN ROWS of the joint fit now
          (`examples/stage13_rows.py` `Rows.lone`, helix up to 4.95 A).
        * FOR THE REVERSED UNIT'S TWO BONDS: `examples/chain_conformers.py frames|settle|collect|polyfind` (lone
          chains of the polymer with every second unit reversed: hh2_tt, hh2_special, hh2_plain on 12 atoms,
          hh4_ff, hh4_hh on 24; DFT in `/opt/dft/chainconf`, `docs/predictions/chain_conformers_pbe_d3bj.json`
          when collected), rows `Rows.reversed_unit`; the candidate's cells on FROZEN geometry, rows
          `Rows.candidate` (`chain_split.py polyfind --freeze FILE`, `crystal_stage.py --frozen FILE`;
          `docs/predictions/chain_split_frozen_belt4.json` is belt4's); six more variables with
          `crystal_stage.py --typed` (`crystalfit.VARIABLES13_TYPED`, 27).
        * THE JOINT FIT WITH THE LINE: DONE, `docs/predictions/crystal_stage_v13_line1.json` (165 -> 95, stopped by
          its own tolerance): lone helix rows to 0.04, gauche alone 50 (DFT 56), repeat alone 4.651 (4.650); gap
          -0.651; alpha b 2.4% short, c 1.3% long; torsion (0.45, -0.52, 0.75), angle stiffness 157.
        * THE REVERSED UNIT'S TWO BONDS: `docs/predictions/chain_stretch_hh_pbe_d3bj.json` (lone head-to-head chain,
          9 DFT points: all-trans is its only minimum; folds held short cost +4.9 to +8.5 and +1.1 to +4.0);
          six typed variables fitted with the rest held (`examples/stage13_subfit.py`, three rounds of fit and
          re-relaxation: `examples/stage13_typed_loop.sh` is only `subfit` then `chain_split.py polyfind --freeze`, repeated):
          `docs/predictions/crystal_stage_v13_typed1.json`, CF2-CF2 (-1.47, 2.01, 0.16), CH2-CH2 (-0.42, 0.60, 1.04).
          "CH2-CH2" is its own lattice torsion type (`fitting.LATTICE_TORSION_TYPES`); every other bond keeps the
          plain triple.
        * PRESET `pvdf-dft-chain1` REGISTERED (`fitting.CRYSTAL13_CHAIN_FITTED_X`, 27 numbers; LATTICE SIDE ONLY, not
          a SimpleFF preset; build polymers inside `with PRESET_PARAMETERS["pvdf-dft-chain1"].applied():` with
          `candidate_screen.make_polymer`, which adds the bond sites; evaluate margins WITH PENDANT TILTS).
          ACCEPTANCE PASSES BOTH HALVES (`docs/predictions/candidate_acceptance_chain1.json`): shear gain -0.275
          (DFT -0.346), all-trans minus TG'TG +0.042 upright / -0.199 relaxed (+0.070 / -0.261), GTG'T minus TG'TG
          +0.242 (+0.197).
        * WHAT chain1 DOES NOT HAVE (5.29's list): alpha's cell (b -2.4%, c +1.3%, gap -0.651); the inside/between
          split of GTG'T minus TG'TG (a CF2-CH2 bond next to a reversed unit is 0.5 too easy to turn: hh_plain row);
          the reversed unit's own backbone angles and bond lengths (its all-trans chain 1.5% long); a fifth of the
          shear gain; absolute cohesion (8 A cutoff); the lone all-trans chain 1.5-3x too stiff when compressed.
        TRAPS: a tilted relaxation (`--dft-gradient`) of a chain with soft turns RUNS AWAY without `--max-move`
        (chain_stretch scan has 0.4 by default; a frame whose predicted energy is far from DFT's, or with a force
        over 0.1 eV/A left, is not a point: `collect --fmax`). A folded conformer that is not a minimum unwinds in
        our potential if its repeat is free: `chain_conformers.py settle --fixed-repeat`. A 36-atom chain in the
        16 A box needs `mlip_frames.py label --pools 1` (memory). `Rows.candidate` needs `--frozen` cells relaxed
        under parameters NEAR the ones being evaluated (first order only): refreeze and refit.
        RUNPOD WORKS AGAIN, THROUGH ITS MCP SERVER (2026-11-18): `https://mcp.runpod.io/` answers from the session
        (the REST host still 502s); the PI saved the key for that host through the secure form, the proxy adds it.
        There are no Runpod tools in the session: speak MCP over HTTP (`initialize`, `tools/list`, `tools/call`; a
        40-line client was in the scratchpad as `runpod/mcp.py`). The key reads pods and capacity (`list-pods`,
        `get-pod`, `get-capacity` = stock and $/h) and creates and deletes pods; it does NOT read the catalog or
        billing (403, final). The old pod `wylsqwj1no08y4` is GONE (404; 0 pods, 0 volumes).
        THE SECOND POD IS GONE AND ITS FIT IS IN THE REPOSITORY (2026-11-19). Pod `dwyejg77awqh50` (NVIDIA A40,
        Secure, $0.59/h) ran 21:40:56 to 00:32:04 UTC: 2.85 h, ABOUT $1.68 (the key cannot read billing; spend so
        far about $17.9 of $46). Removed with `delete-pod` and CONFIRMED BY A READ (`list-pods`: none; `get-pod`:
        404). The queue file never got its `stop` line; nothing reads it now.
        THE FIT `all_s1` (all 558 frames, 485 in the fit, the composition `VAVF|VRVV` held out whole; 240 epochs,
        154 minutes): `docs/predictions/mlip/models/all_s1.model` (14.5 MB, sha256 eabd2345...; saved on the card:
        load with `mlip_second_opinion.file_calculator`, or `mlip_train.py relax --model`), results
        `docs/predictions/mlip_own_all_s1.json`. Held out (69 frames): energy between frames of the same atoms
        1.7 meV per monomer rms (PVDF 0.6 to 1.3, chlorine 3.4, the nitrile 4.5 which is 0.104 kcal/mol and the
        one group outside the 0.1 needed, the composition it never saw 2.9), forces 28.8 meV/A rms (the labels' own rms 1229), stress 0.07 GPa. THE ANCHOR: +0.724 on DFT's
        geometries and +0.725 relaxed in the potential, against DFT's +0.726; volumes of alpha and beta to 0.0%;
        beta's stiffness C_11 23.0 (DFT 22.8), C_22 27.3 (24.5), C_12 7.5 (5.9), C_33 295 (303). The anchor's
        cells are among the frames, so this is the fit reproducing what it was shown, not a prediction.
        THE PRESETS AGAINST DFT ON TrFE AND ON CHLORINE: DONE FOR TWO CELLS (2026-11-19, `docs/MLIP_SECOND_OPINION.md`
        3.10). WHY IT WAS RUN: on the first 14 rows chain1 puts the TrFE matrix 0.17 above belt4 and shrinks the CFE
        pin from +0.26 to +0.04 per fourteen; chain1 had been tested on the reversed-unit candidate only.
        RESULTS (all-trans minus helix, kcal/mol per monomer; cells relaxed in `all_s1`, one DFT single point each,
        0.04 to 0.08 eV/A left): `vdf3-trfe1-4` (VVRV, 25 mol% TrFE) DFT +0.678, chain1 +0.671, belt4 +0.641 (helix
        upright: +0.505, +0.421, +0.577). `vdf2-trfe1-cfe1-4` (VRVC) DFT +0.682, chain1 +0.381, belt4 +1.057: NEITHER
        PRESET HAS CHLORINE. SPLIT BY A LONE CHAIN CUT FROM EACH CELL: VVRV upright +0.505 = +2.935 inside a chain
        - 2.430 between; VRVC +0.682 = +2.979 - 2.298. Between chains both presets have the CFE unit right (0.06,
        0.02); INSIDE THE CHAIN chain1 has it 0.25 too low and belt4 0.28 too high per monomer (about -1.0 and +1.1
        kcal/mol per CFE unit): the bonds at the CFCl carbon take the shared torsion triple. On the TrFE cell belt4
        is right by two errors of 0.7 cancelling; chain1's are +0.13 and -0.22. THE POTENTIAL finds the geometry
        but its own margins are 0.07 to 0.10 from DFT's: relax in it, then ask DFT.
        THE CFCl BONDS HAVE THEIR OWN TORSION TRIPLE NOW: PRESET `pvdf-dft-chain1-cfe` (2026-11-19 05:00 UTC,
        `docs/MLIP_SECOND_OPINION.md` 3.11). Lone chains of the alternating VDF-CFE polymer in DFT
        (`docs/predictions/chain_conformers_cfe_pbe_d3bj.json`): the helix in its two placements -2.475 and -3.189
        kcal/mol per monomer below all-trans; chain1 -1.783 and -2.835. `examples/unit_torsion_fit.py --base
        pvdf-dft-chain1 --type CH2-CFCl --dft THAT --out docs/predictions/unit_torsion_cfe_chain1.json` (25 s) gives
        (+0.4397, -2.1799, +1.5007) and -2.403 / -3.293. Two rows hold ONE combination of the three (the gauche
        well against the trans one); the ridge holds the rest to chain1's triple, so THE BARRIER BETWEEN THE WELLS
        IS NOT FITTED (it matters to the kink barriers, the rate axis). ACCEPTED ON CELLS NOT IN THE FIT
        (`docs/predictions/margins_chain1-cfe_dft_rows.json`): VRVC +0.690 (DFT +0.682; chain1 +0.381), DFT's helix
        placement; VRVK (CTFE) +0.400 (DFT +0.315; chain1 +0.134, belt4 +0.734): inside the 0.1 by 0.015, its
        CF2-CFCl bond has no triple and is trans in both cells; VVRV unchanged +0.671. The CFE cell's split under it:
        +0.715 = +3.161 - 2.445 upright (DFT +0.682 = +2.979 - 2.298), the unit's own part right to 0.05 inside
        and 0.07 between. A MISTAKE CORRECTED: `cfe2_special` and `cfe2_plain` were not conformers (two gauche
        bonds of opposite sense side by side; the chains broke up); dropped, frames in `/opt/dft/cfeconf/broken`.
        Frames `cfe4_tt`, `cfe4_helix`, `cfe4_helix2` (the unit one in four, 24 atoms) are built in
        `/opt/dft/cfeconf` and NOT settled: a further test of the triple, after the unit checks.
        THE NITRILE CELL IN DFT (2026-11-19 06:30 UTC): ALL-TRANS IS BELOW EVERY HELIX PLACEMENT, AND BOTH PRESETS
        HAVE THE SIGN WRONG. VRVA (`vdf2-trfe1-an1-4`): all-trans minus helix -0.386 (TGTG'), -0.404 (TG'TG), -0.678
        (G'TGT); chain1 margin +0.258, belt4 +0.354. Split: DFT -0.386 = +2.128 inside a chain - 2.514 between
        chains; chain1 off by -0.144 inside and +0.590 BETWEEN CHAINS (belt4 -0.277, +0.365). So not a torsion: the
        nitrile's own non-bonded terms. Nitrogen is at UFF (3.66 A, 0.069) in every preset, and the C-N charge
        increment is 0.37 e (a group dipole near 2.1 D; a nitrile's is near 3.9), never a target. A scan
        (scratch): with the increment at 0.75 the three placements go to -0.160, -0.279, -0.547; at 0.99 they
        scatter (cells change basin). `fitting.with_terms(base, lj=, charge_increments=)` and a parameter file's
        `lj` / `charge_increments` put such terms on a preset. ALL THREE PLACEMENTS ARE SPLIT IN DFT
        (`chain_split_dft_vdf2-trfe1-an1-4.json`): TGTG' +2.128 inside -2.514 between; TG'TG +0.212 -0.616; G'TGT
        +1.090 -1.767. chain1 with its cells shear-free (`chain_split.py polyfind --states
        MARGINS_FILE`, which is the right partner of DFT's sheared cells): between +0.81, +0.04, +0.70 off, inside
        -0.16, -0.41, -0.38. LONE CHAINS OF VDF-AN IN DFT (`chain_conformers_an_pbe_d3bj.json`): helix -1.569 and
        +0.367 against all-trans; chain1 -0.535 and +1.286 (both 1.0 high, difference right).
        THE NITRILE IS NOT FIXED AND IS PARKED (2026-11-19 08:30 UTC, docs 3.12): THREE THINGS TRIED, NONE
        ACCEPTED, NO PRESET HAS A NITRILE TERM. (1) C-N increment and nitrogen's wall fitted to the three energies
        (`unit_nonbonded_an_margins_chain1.json`: 0.747 e; 3.70 A, 0.043): energies within 0.05 by two errors of
        0.5 that cancel (+0.50 inside, -0.59 between). (2) The same fitted to the between-chains parts
        (`unit_nonbonded_fit.py --between`, `unit_nonbonded_an_between_chain1.json`: 0.632 e; 3.445 A, 0.052): the
        parts stay 0.30 and 0.17 off (the cells change basin as the charge moves: 0.19 for 0.002 e), then a
        torsion triple for CH2-CHCN on the lone chains (`unit_torsion_an_chain1.json`, coefficients of 3). TESTS:
        VRVA's margin -0.185 without the torsion and +0.046 with it (DFT -0.386); THE SECOND CELL VVVA
        (`vdf3-an1-4`, DFT -0.091 = +0.878 inside - 0.969 between) -0.549 and +0.164. (3) An N-H contact term
        (`with_terms(lj_pairs=)`, a scratch scan): moves the between-chains parts by 0.15 at most. chain1 on VVVA:
        -0.292 on another placement, -0.751 on DFT's; belt4 +0.058. chain1's errors on the two cells do not agree
        (VVVA -0.53 inside, -0.15 between; VRVA -0.16..-0.41 inside, +0.04..+0.81 between), so there is no one
        correction of this form. WHAT IT NEEDS (a guess, not tried): a charge model for the nitrile that moves
        charge on the carbon it sits on and that carbon's hydrogen, refitted on nitrile cells and lone chains.
        THE POTENTIAL HAS BOTH NITRILE CELLS TO 0.03 AND 0.08: for the 7 nitrile rows it is the evaluator.
        IF THE UNIT PINS ALL-TRANS AT ONE IN FOURTEEN AS IT DOES AT ONE IN FOUR (-3.3 to -4.3 kcal/mol per unit),
        THOSE ROWS MOVE DOWN TOWARD THE WINDOW, WHERE THE CFE ROWS MOVE UP (terpolymer +0.359 under chain1-cfe,
        matrix +0.222). THAT IS THE QUESTION THE POD RUN PUT TO THE PI WOULD ANSWER (all 35 rows' cells relaxed in
        the potential; no answer yet, NO POD WITHOUT IT).
        FOUR FOURTEEN-MONOMER ROWS IN THE POTENTIAL (2026-11-19 13:55 UTC, docs 3.14,
        `docs/predictions/rows14_potential.json`; NO DFT AT THIS SIZE): matrix `vdf10-trfe4-14` -0.024 (chain1-cfe
        +0.222, belt4 +0.056); CFE terpolymer `vdf9-trfe4-cfe1-14` +0.078 (+0.359, +0.317); nitrile alone
        `vdf9-trfe4-an1-14` -0.118 (belt4 +0.197); THE belt4 LEADER `vdf9-trfe3-an1-rvdf1-14` -0.055 (belt4
        +0.037): at the foot of the window, NOT CONFIRMED AND NOT THROWN OUT (the caution at the head of
        `docs/CANDIDATE_MATERIALS.md` says so). The potential and chain1-cfe AGREE on what a CFE does (+0.10,
        +0.14) and are 0.25 APART ON THE MATRIX'S LEVEL. The nitrile at one in fourteen: -1.3 kcal/mol per unit
        (all-trans pin) in the potential; belt4 +2.0. THE POTENTIAL'S OWN MARGIN: within 0.11 of DFT on six small
        cells, 0.36 OFF ON THE TFE CELL. How: cells from `rows_under_preset.py` (nitrile rows:
        `/opt/dft/cells_chain1-cfe_an.json`, geometries only) written by `row_cells_frames.py write --kmesh "4 2
        1"` into `/opt/dft/rows14an` and `/opt/dft/rows14`; `mlip_train.py relax --model
        docs/predictions/mlip/models/all_s1.model --name all_s1 --cells ... --threads 4 --fmax 0.02 --steps 800`
        (five to ten minutes a cell); `row_cells_frames.py potential --frames DIRS --rows ROWS --out JSON`.
        WHICH LEVEL IS RIGHT: ANSWERED BY DFT (2026-11-19 15:40 UTC, docs 3.15). `vdf4-trfe2-6`, VVRVVR, the
        matrix's motif, 72 atoms, its four helix placements one cell by symmetry: DFT -0.044 (= +1.603 inside
        a chain -1.647 between chains; 0.09 eV/A left), THE POTENTIAL -0.096, belt4 -0.014 (by two errors of 0.79
        and 0.76 that cancel), chain1 +0.139 shear free (upright +0.078 = +2.004 -1.927). DFT IS WITH THE
        POTENTIAL; chain1 IS 0.18 HIGH. The gap chain1 minus potential is 0.235 here and 0.246 on the
        fourteen-monomer matrix. chain1's TWO PARTS ARE OFF ON EVERY CELL, the inside of the chain +0.13 to
        +0.40 and between chains -0.15 to -0.29; its passes on the four-monomer cells are where the two are
        equal (VVRV +0.13/-0.22, VRVC +0.18/-0.15, VRVK +0.30/-0.15, VRVF +0.26/-0.29; VVRVVR +0.40/-0.28).
        WHAT IS NEW IN VVRVVR: each chain holds one TrFE on each gauche sense (DFT's helix: CF2-CHF bonds at +65
        and -56 degrees); VVRV's lowest helix, the only one DFT saw, has the favourable sense only. A torsion
        term is even in the angle and cannot tell the senses apart (BUT SEE THE LONE CHAINS BELOW: the error is
        not on the other sense). SO THE 25 ROWS' LEVELS UNDER chain1-cfe READ UP TO 0.2 HIGH at 29 mol% TrFE and
        above (right at one TrFE in four, 0.18 high at one in three); the differences a unit makes stand.
        ESTIMATES (one cell's offsets carried to fourteen monomers, NO DFT there): matrix `vdf10-trfe4-14` +0.03
        (from the potential +0.05), +0.04 (from chain1-cfe -0.18), +0.03 (from belt4 -0.03): just under the
        window; CFE terpolymer `vdf9-trfe4-cfe1-14` +0.13 and +0.18: at the window's top. A FOURTH (chain1-cfe
        with the TrFE term: matrix +0.100, terpolymer +0.239, `margins_trfe_torsion_rows14.json`; corrected by
        its errors on the two TrFE cells +0.08 to +0.20) WIDENS IT: THE MATRIX IS SOMEWHERE FROM +0.03 TO +0.2,
        the terpolymer 0.10 to 0.14 above it. Only DFT on those cells narrows that. Command, for the
        record: `python examples/unit_dft_check.py --row vdf4-trfe2-6 --work /opt/dft/unit4 --model
        docs/predictions/mlip/models/all_s1.model --pw /opt/mm/qe/bin/pw.x --pseudo /opt/pseudo --check
        /opt/dft/chainsplit/check --threads 4 --np 4` (1.9 hours). NOTE: "DFT" on a row's cell is a DFT single
        point on the cell the potential relaxed (forces left are stored), as in docs 3.6 and 3.10.
        THE FOUR AXES AT THE POTENTIAL'S MARGINS (2026-11-19 15:50 UTC; `examples/four_axes_at_margins.py --out
        docs/predictions/four_axes_potential_margins.json`; table at the head of `docs/CANDIDATE_MATERIALS.md`;
        test `tests/test_four_axes_at_margins.py`): the potential's margins of the four fourteen-monomer rows put
        through `pin_design.py`, EVERYTHING ELSE belt4's (dipole, switch strain, KINK ENERGIES, stiffness: no DFT
        test on any). (1) THE CELL OF THE FILM'S COMPOSITION LANDS ON THE MEASURED FILM: at +0.078, tan delta
        0.088 and AC strain at the heat budget 0.081% against the film's calibrated 0.060 and 0.066% (0.03 to
        0.09 and 0.055 to 0.081% over margin -0.1 to +0.1); at belt4's +0.317 it had 0.004%. (2) THE LOSS MODEL
        TAKES A MARGIN BELOW ZERO AS ZERO; below zero the cell is a ferroelectric the model does not hold. (3) THE
        LEADER (-0.055 in the potential): if not below zero it leads on all four axes as before (2.6 to 6.9 MHz,
        tan delta 0.003 to 0.012, 0.133 to 0.168%, 7.2 to 8.8 pN/nm2 through the crystal against the film's 364
        kHz, 0.060, 0.066%, 6.5); if below, it is not a candidate of this kind. ITS STANDING IS ONE QUESTION: THE
        SIGN OF ITS MARGIN, which the potential's 0.1 cannot answer. WHAT SETTLES IT: DFT on its all-trans and
        TGTG' cells as the potential relaxed them (170 atoms each,
        `/opt/dft/rows14an/row_vdf9-trfe3-an1-rvdf1-14_{TTTT,TGTGp}_cfe_in_all_s1.start.xyz`; an estimated one to
        two hours a cell on four cores, 15 GB of memory here: NOT STARTED). (4) The nitrile alone (-0.118) is
        below zero over all of the bar but its last 0.02. (5) THE KINK ENERGIES ARE THE WEAKEST INPUT: belt4's,
        untested, and belt4 has the nitrile wrong in sign on the margin.
        LONE VDF-TrFE CHAINS, DONE (16:12 UTC, docs 3.16, `chain_conformers_trfe_pbe_d3bj.json`; `trfe2` in
        `chain_conformers.CHAINS`; `/opt/dft/trfeconf`): helix against all-trans per monomer, DFT -2.802 (CF2-CHF
        at -62, favourable) and +0.631 (+49); chain1 -3.071 and +0.551; belt4 -3.741 and +0.439. chain1 HAS THE
        DIFFERENCE BETWEEN THE SENSES (3.62 against 3.43); ITS ERROR IS ON THE FAVOURABLE ONE (0.27 too low), so
        the suspicion that the other sense carried VVRVVR's error IS WRONG. A CF2-CHF TRIPLE FITTED TO THE TWO ROWS
        (`unit_torsion_fit.py --base pvdf-dft-chain1-cfe --type CF2-CHF`, (+0.079, +0.400, +0.326),
        `docs/predictions/unit_torsion_trfe_chain1.json`) AND TRIED ON FIVE CELLS (`rows_under_preset.py --preset
        THAT_FILE --rows ... --polished docs/predictions/polished_margins_belt4_dft_rows.json --references
        docs/predictions/dilute_kink_belt4_dft_rows.json`; `margins_trfe_torsion_dft_rows.json`): off against DFT
        VVRV -0.10 (was -0.01), VVRVVR +0.02 (was +0.18), VRVC -0.09 (+0.01), VRVK -0.05 (+0.09), VRVF -0.13
        (+0.02); root mean square 0.09 with it and 0.09 without. NOT ADOPTED, NO PRESET: with the chain right the
        error between chains (all-trans 0.15 to 0.29 too strongly bound) is no longer cancelled, and that one is
        in the walls and charges every unit shares: a refit of the preset, not a unit's term. SCORE ON THOSE FIVE
        CELLS, largest and root mean square: the potential 0.36 (TFE) and 0.18; belt4 0.42 and 0.29; chain1-cfe
        0.18 and 0.09. NO EVALUATOR IS WITHIN 0.1 ON ALL FIVE.
        RUNNING (16:27 UTC, one core, ends by itself): the fourteen-monomer matrix and CFE terpolymer under the
        fitted TrFE term, for a third estimate of their level: `rows_under_preset.py --preset
        /opt/dft/unit_torsion_trfe_chain1.json --redo-failed --minutes 100 --rows
        vdf10-trfe4-14,vdf9-trfe4-cfe1-14 --out /opt/dft/margins_trfe_torsion.json`.
        A KINK'S ENERGY AGAINST DFT, DONE FOR ONE CELL (18:10 UTC, docs 3.17,
        `docs/predictions/kink_dft_vdf4-trfe2-6.json`, 4 labelled frames in
        `mlip/kink_vdf4-trfe2-6_frames_pbe_d3bj.extxyz`; the kink energies set frequency and tan delta and had no
        test). The six-monomer matrix cell, both chains kinked (T G+ T G- T at bonds 0 to 2), lattice held; a kink
        pair does not fit a six-monomer repeat (gauche bonds flatten to about 150 and -140; IT COLLAPSES WHEN THE
        CELL IS FREED), so this is a squeezed kink, the same structure in every evaluator. kcal/mol per kink pair:
        AT belt4's ALL-TRANS LATTICE (53.5 A3 a monomer, 7% too dense) DFT 6.01, the potential 6.68, belt4 6.60;
        AT THE LATTICE THE POTENTIAL RELAXES THE ALL-TRANS CELL TO (57.5 A3) DFT 4.15, the potential 4.75; chain1
        at its own lattice (57.8 A3) 3.66 (the potential there 5.01, no DFT). SO: belt4's KINK IS 59% DEAR against
        DFT at DFT's lattice (10% at its own); chain1 12% CHEAP; THE POTENTIAL 0.60 TO 0.67 HIGH AT BOTH (never
        shown a kink) with the lattice's part right (1.92 against 1.85). NOT TESTED: the ratios between cells,
        which are what the rates are. How: `python examples/dilute_kink.py --row vdf4-trfe2-6 --nx 1 --ny 1
        --kinked all --stages flexible_host,c_free,cell_free [--preset pvdf-dft-chain1] --out
        /opt/dft/kink6/dilute_kink_{belt4,chain1}.json`; `python examples/row_cells_frames.py kink --dilute THAT
        --key "vdf4-trfe2-6@0 1x1 kinked 0,1" --stages flexible_host --kmesh "4 2 2" [--tag chain1] --out
        /opt/dft/kink6`; `mlip_train.py relax ... --mask "0 0 0 0 0 0" --fmax 0.02` (atoms only);
        `mlip_frames.py label --frames /opt/dft/kink6 --names ... --np 3` (25 minutes a cell);
        `row_cells_frames.py kink-collect --pairs NAME=TT.dft.xyz,KINKED.dft.xyz ... --out JSON --extxyz FILE`.
        KINK SHEETS OF THE FOUR FOURTEEN-MONOMER CELLS IN THE POTENTIAL: SEVEN OF 112 SITES DONE (20:00 UTC, docs
        3.18, `docs/predictions/sheets14_potential.json` and `sheets14_potential_more_sites.json`; cell free,
        each from belt4's cell for that site, every kink kept; NO DFT on these cells; the potential is 0.6 high
        on the one kink DFT has seen). kcal/mol per kink pair, the potential (belt4): matrix site 12~ 4.54 (4.66),
        site 12 3.59 (4.81); CFE cell site 10~ 4.59 (3.42), site 0 4.16 (3.95), site 8~ 6.60 (3.57); nitrile alone
        site 8~ 2.95 (2.15); leader site 8~ 2.49 (1.39). SITE BY SITE THE TWO DIFFER BY -1.2 TO +3.0: belt4's
        PROFILE ALONG THE SEQUENCE IS NOT CONFIRMED (a first reading on four cells, "leader -2.10 against belt4's
        -2.02, speed confirmed", was of one site each and WAS WITHDRAWN). WHAT STANDS: the order of the lowest
        sheets found, leader 2.49 < nitrile alone 2.95 < matrix 3.59 < CFE cell 4.16; THE CFE's CUT OF THE KINK
        SHEET (1.24 in belt4) IS NOT THERE. THE FOUR AXES WITH EACH CELL'S SITES FLOORED AT THE LOWEST FOUND
        (`four_axes_at_margins.py --sheet-shift vdf10-trfe4-14=-1.067 vdf9-trfe4-cfe1-14=0.744
        vdf9-trfe4-an1-14=0.803 vdf9-trfe3-an1-rvdf1-14=1.092 --out
        docs/predictions/four_axes_potential_margins_sheets.json`; second table at the head of
        `docs/CANDIDATE_MATERIALS.md`): leader 3.8 MHz, tan delta 0.006, 0.145%, 7.9 pN/nm2 through the crystal;
        nitrile alone 1.6 MHz, 0.013, 0.123%, 8.5; matrix 697 kHz, 0.027, 0.086%, 9.8 (each IF ITS MARGIN IS NOT
        BELOW ZERO); film as calibrated 364 kHz, 0.060, 0.066%, 6.5. PROVISIONAL: THE RATES NEED EVERY SITE IN
        THE POTENTIAL (14 sites, two hands, four cells; a quarter of an hour each here: A JOB FOR THE POD). How:
        `python examples/row_cells_frames.py kink --scan docs/predictions/kink_site_scan_belt4.json --key
        "ROW@0~" --sites N --kmesh "4 2 1" --out DIR`; `mlip_train.py relax --model
        docs/predictions/mlip/models/all_s1.model --name all_s1 --cells DIR/sheet_*.start.xyz --out DIR --threads
        3 --fmax 0.02 --steps 800`; `row_cells_frames.py sheets-potential --frames DIR --out JSON` (cells in
        `/opt/dft/sheets14` and `/opt/dft/sheets14b`).
        MATRIX AND TERPOLYMER UNDER THE TrFE TERM (not adopted): +0.100 and +0.239
        (`docs/predictions/margins_trfe_torsion_rows14.json`; chain1-cfe +0.222 and +0.359, the potential -0.024
        and +0.078).
        ETHYLENE (16:50 UTC): VRVE DFT -0.146 (= +1.007 -1.153), chain1 -0.232 (PASSES AT THE EDGE; parts +0.04,
        -0.23), belt4 +0.500 (WRONG SIGN), the potential -0.082. VF (20:37 UTC): VRVY DFT +0.290 (= +2.081
        -1.791; DFT LEAVES 0.32 AND 0.35 eV/A on the two cells, the most of any pair), chain1 +0.479 (FAILS BY
        0.19; parts +0.13, +0.07), belt4 +0.452, the potential +0.274. DFE (21:08 UTC): VRVX DFT +0.961 (=
        +2.213 -1.251; 0.29 eV/A left on both cells), chain1 +0.955 (PASSES CLEANLY; parts +0.04, -0.07), belt4
        +1.234, THE POTENTIAL +0.309: 0.65 LOW, ITS WORST CELL. THE UNIT CHECKS ARE COMPLETE: every unit of the 35
        rows has a DFT cell. THE POTENTIAL MISSES TFE (0.36) AND DFE (0.65) and is within 0.11 on the other ten
        cells; chain1-cfe is within 0.09 on seven, 0.18 to 0.19 off on VVRVVR and VF, right by cancellation on
        CDFE, wrong on the nitrile. USE EACH WHERE IT HAS PASSED. Rows by evaluator: nitrile rows (7) the
        potential; TFE (2) and DFE (1) rows chain1-cfe; VF (1) the potential; the rest either, with chain1-cfe's
        levels up to 0.2 high at 29 mol% TrFE. TODAY'S LABELLED CELLS (the `row_cells_*_frames_pbe_d3bj.extxyz`
        of the eight small rows, the lone chains, `kink_vdf4-trfe2-6_frames_pbe_d3bj.extxyz`) ARE NOT IN THE
        POTENTIAL'S TRAINING SET: a refit with them (2.85 hours on the A40 last time) comes before the potential
        is used on TFE or DFE rows or on kinks.
        THE PRESET'S TABLE IS NOW 27 ROWS (22:05 UTC, `docs/predictions/margins_chain1-cfe.json`): added
        `vdf9-trfe4-e1-14` -0.046 (ethylene -3.8 kcal/mol per unit: a pin for all-trans) and `vdf9-trfe4-dfe1-14`
        +0.287 (DFE +0.9 per unit). Not in it: the seven nitrile rows and the VF row.
        FOUR MORE KINK-SHEET SITES IN THE POTENTIAL (22:31 UTC, `docs/predictions/sheets14_potential_third_batch.json`,
        cells in `/opt/dft/sheets14c`; ELEVEN OF 112 NOW): leader `@0` site 18 2.41 (belt4 1.86), LEADER `@0~`
        SITE 18 0.65 (belt4 1.94; kink kept, gauche -72 and +61), nitrile alone `@0` site 0 3.40 (3.25), NITRILE
        ALONE `@0~` SITE 10: THE KINK DOES NOT HOLD (back to all-trans; belt4 has a sheet there at 2.26). LOWEST
        FOUND NOW: leader 0.65, nitrile alone 2.95, matrix 3.59, CFE cell 4.16. A sheet for 0.65 is a kink pair
        for kT: outside the rate model's dilute range (the formula gives 83 MHz; on the leader's next site, 2.41,
        4.4 MHz, tan delta 0.005, 0.149%, 8.1 pN/nm2 through the crystal). `four_axes_potential_margins_sheets.json`
        is recomputed with the leader's shift -0.740 (= 0.655 - 1.395); this SUPERSEDES the 3.8 MHz row written
        at 20:00 UTC above. ONE RELAXATION FROM ONE STARTING CELL, NO DFT: the full scan and DFT on the lowest
        sheets are the pod's job.
        A REAL KINK SHEET TO DFT, RUNNING (22:42 UTC, four cores, one pool, 6.6 GB by pw.x's estimate, ends by
        itself): PVDF, ten monomers, 120 atoms, the kinked crystal in its own cell (the six-monomer cell holds
        only one squeezed kink: its other five sites and the other hand all collapse under chain1). belt4 +5.298
        kcal/mol per kink pair (5.063 with tilts); THE POTENTIAL +4.685, both cells relaxed with the cell free,
        kink kept (81, 178, -68). `python examples/row_cells_frames.py kink --sheet-of "vdf10@0" --dilute
        docs/predictions/dilute_kink_belt4.json --kmesh "4 2 1" --out /opt/dft/kink10` (new option: a stored
        kinked crystal of `kink_migration_belt4.json` and its all-trans cell); `mlip_train.py relax ... --cells
        /opt/dft/kink10/sheet_vdf10_belt4_{TT,ref}.start.xyz --out /opt/dft/kink10 --threads 4 --fmax 0.02`;
        `python -I examples/mlip_frames.py label --frames /opt/dft/kink10 --names
        sheet_vdf10_belt4_ref_in_all_s1,sheet_vdf10_belt4_TT_in_all_s1 --pw /opt/mm/qe/bin/pw.x --pseudo
        /opt/pseudo --np 4 --pools 1` (the same command goes on after a stop). DONE 23:49 UTC (33 minutes a
        cell): DFT'S SINGLE POINTS ON THE TWO CELLS GIVE +4.820 (the potential +4.685 on the same cells, belt4
        +5.298 in its own). BUT DFT LEAVES 0.16 eV/A AND 0.19 GPa ON THE KINKED CELL (0.01 and 0.08 on the
        all-trans one), and one Newton step toward DFT's own minimum (`mlip_train.py relax --dft-gradient
        --max-move 0.5 --cells /opt/dft/kink10/sheet_vdf10_belt4_{TT,ref}_in_all_s1.dft.xyz --out
        /opt/dft/kink10/tilted`) moves atoms up to 0.57 A, shears the cell by 0.8 degrees and EXPECTS DFT AT
        +4.03. That is too large a step to quote unchecked: DFT ON THE STEPPED CELL IS RUNNING (23:56 UTC, about
        35 minutes; `mlip_frames.py label --frames /opt/dft/kink10/tilted --names
        sheet_vdf10_belt4_ref_in_all_s1_in_all_s1 --pw /opt/mm/qe/bin/pw.x --pseudo /opt/pseudo --np 4 --pools 1`,
        the same command goes on after a stop). THEN `row_cells_frames.py kink-collect --pairs
        vdf10=/opt/dft/kink10/sheet_vdf10_belt4_TT_in_all_s1.dft.xyz,/opt/dft/kink10/sheet_vdf10_belt4_ref_in_all_s1.dft.xyz
        --steps vdf10=/opt/dft/kink10/tilted/sheet_vdf10_belt4_ref_in_all_s1_in_all_s1.dft.xyz --also
        belt4_with_tilts=5.063 --out docs/predictions/kink_dft_vdf10.json --extxyz
        docs/predictions/mlip/kink_vdf10_frames_pbe_d3bj.extxyz` (`--steps` is new: the kinked cell after each
        Newton step, labelled), docs 3.19, a test. The same step on the six-monomer kink of 3.17 (lattice
        held) expects DFT to drop by only 0.04 to 0.05: those numbers stand as they are.
        THE STEPPED CELL'S DFT IS IN (00:34 UTC): +4.143 (expected +4.027), 0.06 eV/A and 0.18 GPa left, the
        kink kept (80, -179, -69). Collected with `--steps` into `docs/predictions/kink_dft_vdf10.json`; test
        `test_a_real_kink_sheet_against_dft_needs_dfts_own_geometry`. THE CAUSE IS THE STOP, NOT ONLY THE
        POTENTIAL (01:10 UTC): at `--fmax 0.02` a cell of 120 atoms does not answer a stress under 0.35 GPa
        (the cell's share of the stop is the stress times the volume over the atom count), and the potential's
        own stress on its "relaxed" kinked cell is 0.29 GPa (0.07 on the all-trans one). On DFT's stepped
        cell the plain potential has the sheet at +4.431, BELOW its own +4.685: the cell it stopped at was not
        its minimum. Same geometry, like for like: the potential 0.14 low on its loose cell, 0.29 high on the
        stepped one. EVERY FOURTEEN-MONOMER SHEET AND ROW OF 3.14 AND 3.18 WAS STOPPED AT 0.02 TOO: upper
        bounds by an amount not yet sized. A second tilted step at `--fmax 0.003` runs away along a soft
        direction (a +8.4%): discarded (`/opt/dft/kink10/tilted2`). RUNNING (01:11 UTC, about 20 minutes a
        cell): the plain potential to `--fmax 0.003 --steps 420` on both cells, `--out /opt/dft/kink10/tight
        --frame sheet_vdf10_belt4_{ref,TT}_tight` (the kinked one from the stepped cell,
        `tilted/sheet_vdf10_belt4_ref_in_all_s1_in_all_s1.start.xyz`). THEN DFT on the tight kinked cell (35
        minutes; and on the all-trans one if it moved), `kink-collect` again, docs 3.19. THE CONTAINER WAS
        RESTARTED AT 01:10 UTC: files kept, running jobs lost (only relaxations in the potential).
        THE TIGHT STOP FINDS A DIFFERENT SHEET (01:29 UTC): in the plain potential the kinked cell keeps
        shearing along the chains, beta 91 -> 111.4 and alpha 90.4 -> 97.9 degrees, a +10%, volume -0.9%,
        the kink kept and sharper (62, 172, -59), and the sheet's energy falls from +4.685 to +0.846 kcal/mol
        per kink pair, NOT CONVERGED after 420 steps (`/opt/dft/kink10/tight/sheet_vdf10_belt4_ref_tight`).
        The all-trans cell does not move (29 steps, -0.1 meV per monomer). The shear is the chains next to
        each other along a sliding along their axes, about 2.0 A, most of a monomer (2.49 A): the kinks of
        neighbouring chains are no longer at one height, THE SHEET IS INCLINED. DFT's stress on both cells it
        has seen points the same way, undiminished after the first step (xz +0.19, +0.18 GPa: about 0.25
        kcal/mol per kink pair per degree of beta; the potential's own slope there is 0.40). BELT4'S KINK
        SHEETS ARE ALL UPRIGHT (its scans start at alpha = beta = 90 and its kink cells have no shear along
        the chains), AND TEN OF THE ELEVEN FOURTEEN-MONOMER SHEETS OF 3.18 STALLED NEAR UPRIGHT (c tilted
        under 1.2 A). THE ELEVENTH IS THE LEADER'S 0.65: its cell had sheared 3.1 A along a and 4.0 A along b
        (399 steps). SO 0.65 IS NOT THE LEADER'S PROPERTY, IT IS THE ONE RELAXATION THAT GOT PART WAY DOWN A
        SLOPE EVERY SHEET IS LIKELY TO HAVE; the order of the lowest sheets in 3.18 and the frequencies made
        from it (leader 83 MHz / 4.4 MHz, nitrile 1.6 MHz, matrix 697 kHz) ARE WITHDRAWN until the sheets are
        relaxed tight in every cell. RUNNING (01:31 UTC): the kinked cell carried on, `--cells
        /opt/dft/kink10/tight/sheet_vdf10_belt4_ref_tight.start.xyz --fmax 0.003 --steps 480 --frame
        sheet_vdf10_belt4_ref_tight2`. THEN DFT on that cell (35 minutes): is the drop DFT's too.
        THE INCLINED SHEET IS BELOW ALL-TRANS, IN EVERY EVALUATOR (02:50 UTC). The potential converged
        (480 more steps): cell [5.42, 8.87, 24.83, 102.5, 113.6, 83.1], which is a - c0 and b - c0 of the
        all-trans lattice (each neighbour's kink one monomer along), 53.3 A3 a monomer, kink 60, 172, -58, the
        rest of the chain exactly trans: -1.63 kcal/mol per kink pair against the all-trans cell (tight;
        all-trans itself is at its minimum upright: a control from eight slides in the primitive cell,
        `/opt/dft/kink10/tt_shear_potential.json`, finds nothing lower and returns the same crystal from the
        one-monomer slide). DFT ON THAT CELL: SCF energy -2758.08253544 Ry at 6e-8 Ry accuracy against the
        all-trans cell's -2758.07005033: -1.96 KCAL/MOL PER KINK PAIR. The run was cut by its own time limit
        (60 minutes, it shared the cores) before forces: RE-RUNNING on three cores for the label (02:52 UTC,
        about 50 minutes; `mlip_frames.py label --frames /opt/dft/kink10/tight --names
        sheet_vdf10_belt4_ref_tight2 --np 3 --pools 1`). POLYFIND ITSELF, WITH THE SHEAR IT ALREADY HAS
        (`examples/kink_sheet_shear.py`, new: a stored kink sheet in `MixedModel(shear=True)`, from the upright
        cell and from starts a monomer along a and b): belt4 upright +5.298 (a local minimum: a kick returns
        +5.276), INCLINED -1.172 at shear (-0.86, -0.84) monomer; other staggers +0.4 to +4.4; chain1 -0.781.
        EVERY KINK CELL OF dilute_kink.py, kink_migration.py AND kink_site_scan.py IS UPRIGHT
        (`MixedModel` without `shear=True`): the frequency and loss axes of every candidate are made of
        upright sheets. Consistent with PVDF's margin: its helix is 0.73 kcal/mol per monomer below all-trans
        and a kink pair is one turn of that helix, two monomers. RUNNING (02:52 UTC, one core): the same for
        the leader's and the matrix's stored sites, `--key "vdf9-trfe3-an1-rvdf1-14@0~" --sites 18,8` and so
        on into `/opt/dft/kink10/shear_rows_belt4.json` (to be copied to
        `docs/predictions/kink_sheet_shear_belt4.json` when done; log `/opt/dft/kink10/kink_sheet_shear_rows.log`);
        the CFE and nitrile cells' sites when DFT frees the cores. THE 112-SITE SCAN WITH THE SHEAR FREE IS A
        JOB FOR THIS MACHINE IN POLYFIND (a minute or two a start), not for a pod.
        DFT'S LABEL OF THE INCLINED CELL IS IN (03:35 UTC): -1.959 kcal/mol per kink pair (the potential
        -1.655 on the same pair), 0.07 eV/A and 0.04 GPa left: DFT is at rest on the potential's cell.
        Stored: `docs/predictions/kink_dft_vdf10.json` (pairs `vdf10`, `vdf10_inclined`; `kink-collect
        --moved vdf10_inclined` leaves polyfind's upright number off it), 4 frames in the extxyz, tests
        `test_a_real_kink_sheet_against_dft_needs_dfts_own_geometry` and
        `test_the_inclined_kink_sheet_is_below_all_trans_in_dft`; docs 3.19, notices in 3.17 and 3.18, the
        head of `docs/CANDIDATE_MATERIALS.md` (caution of 2026-11-20; the table with the potential's sheets
        marked withdrawn). PVDF's belt4 runs: `docs/predictions/kink_sheet_shear_belt4.json`.
        FIRST CANDIDATE SITES WITH THE SHEAR FREE (belt4, nine of 112): leader `@0~` 18: 1.935 -> 1.991,
        `@0~` 8: 1.395 -> 1.358, `@0` 18: 1.856 -> 1.895 (NONE INCLINES: every start returns upright); matrix
        `@0` 12: 4.815 -> 4.091, `@0~` 12: 4.659 -> 1.136; film's cell `@0` 0: 3.952 -> -0.131, 2: 5.297 ->
        3.450, 26: 4.509 -> 1.124, 4: 4.559 -> 1.099. THE ORDER OF THE CELLS REVERSES.
        THE SCAN IS DONE (06:01 UTC 2026-11-20): ALL 108 VALID SITES OF THE FOUR CELLS, belt4, shear free, the
        upright start and four starts a monomer along a and b. Stored: `docs/predictions/
        kink_sheet_shear_belt4_sites.json` (every site's starts with states; `kink_sheet_shear.py --merge
        /opt/dft/kink10/shear_scan_*.json --gather ...`), `kink_site_scan_belt4_shear.json` (`--as-scan`: the
        scan with each site at its lowest sheet), `four_axes_potential_margins_shear.json`
        (`four_axes_at_margins.py --sites docs/predictions/kink_site_scan_belt4_shear.json --out ...`); tests
        in `tests/test_four_axes_at_margins.py`; docs 3.20 and the head of `docs/CANDIDATE_MATERIALS.md`.
        CHEAPEST SHEET PER CELL, upright -> free, kcal/mol per kink pair: film cell +3.42 -> -2.27 (at V-CFE,
        shear (0, -1.8) monomers; 9 of 28 sites below all-trans); nitrile cell +2.15 -> -2.08 (6 of 26);
        LEADER +1.40 -> +0.65 (0 of 26; 20 of its 26 sites stay upright); matrix +4.66 -> +0.94 (0 of 28).
        THE ORDER REVERSES. FOUR AXES AT THE POTENTIAL'S MARGINS ON THESE SHEETS (34 kHz, heat budget;
        frequency by the sum of exp(-E/kT) and by the sites' share 1/(1+exp(E/kT)), both in
        `pin_design.py` now): film cell 364 kHz, tan delta 0.088, 0.081%, 7.4 pN/nm2; nitrile alone 306 to
        335 kHz, 0.060 to 0.066, 0.066 to 0.068%, 4.5 to 4.7; LEADER 4 to 30 kHz, 0.65 to 0.85, 0.007 to
        0.020%, 0.4 to 1.1; matrix 5 to 34 kHz, 0.53 to 1.12, 0.008 to 0.021%, 0.9 to 2.4. THE LEADER IS WITH
        THE MATRIX ON THESE SHEETS AND NONE OF THE FOUR BEATS THE FILM. NOT A RESULT ON THE MATERIALS YET.
        WHY THE LEADER'S STAY UPRIGHT: the shear by two monomers along b that the others' sheets take costs
        its ALL-TRANS cell 7.4 kcal/mol per cell (film cell 1.7, matrix 2.3, nitrile cell 2.5); its sheets at
        that shear are 2.5 to 3.9 above the nitrile cell's at like sites, the same difference per chain. The
        upright all-trans cell is the lowest in all four from twelve whole-monomer starts (film 1.7 to 5.3
        above, matrix 2.3 to 8.2, leader 4.4 to 11.7; scratch `tt_registry.py` and `shear_polish.py --starts`
        into `/opt/dft/kink10/tt_registry_*.json`, `registry_margins_ROW.json`). Pendant tilts change little
        (film -2.27 -> -2.36; leader 0.65 -> 0.50; leader's sheared all-trans 7.4 -> 7.6). BELT4 IS WEAKEST
        HERE (3.9: three times DFT's cost for sliding chains; wrong sign on the reversed-unit candidate's
        sheared all-trans cell). chain1-cfe: film cell's all-trans at (0, +-2) +0.9 per cell, leader's (0, +2)
        +6.1, most other starts slide back upright. A PERIODIC CELL PUTS EVERY UNIT BESIDE ITS LIKE; a
        statistical copolymer cannot: what the leader's all-trans costs out of register is open.
        THE MARGINS DO NOT MOVE (07:30 UTC): `shear_polish.py --rows ROW --starts "0:-2,0:2,-1:-1,1:1,1:-1,
        -1:1,2:0,-2:0,0:-1,0:1,1:0,-1:0" --out /opt/dft/kink10/registry_margins_ROW.json` for the film cell
        and the leader, every cell from twelve whole-monomer shears: rigid margins with the shear free +0.2511
        and +0.1627, as `sheared_margins_belt4.json` has them; one helix conformer of the leader's that does
        not set the margin (TG'TG) found a cell 1 kcal/mol lower. Matrix and nitrile cell: all-trans only.
        chain1-cfe ON THE KINK SITES: its own sheared starts SLIDE BACK UPRIGHT (its sliding cost is a third
        of belt4's), so `kink_sheet_shear.py --preset pvdf-dft-chain1-cfe --starts "" --from-sheets
        /opt/dft/kink10/belt4_sites_copy.json` starts each site from the sheets belt4 found. Film cell `@0`:
        site 12 (V-CFE) -1.65 (belt4 -2.27), 10: -0.94 (-1.12), 0: -0.53 (-0.13): ITS SHEETS ARE BELOW
        ALL-TRANS THERE TOO. Leader `@0` site 10: +1.66 (belt4 +0.65), its upright sheet loses the kink and two
        stored sheets collapse (chain1-cfe has no nitrile term): stopped after that site. RUNNING (job
        `bb0ps4aqk`, since 07:00): the film cell's sites "12,10,0,4,16,6,22" of `@0` and "2,14,12,0" of `@0~`
        into `/opt/dft/kink10/shear_from_belt4_chain1cfe_vdf9-trfe4-cfe1-14.json`.
        chain1-cfe, film cell, nine sites (stored `docs/predictions/kink_sheet_shear_chain1cfe_sites.json`):
        EIGHT BELOW ALL-TRANS, -0.36 to -1.65 and -4.50 at `@0~` site 14 (CFE-V; belt4 -1.17).
        THE POTENTIAL AT A LOOSE STOP IS NO USE HERE: the leader's upright all-trans cell relaxed to 0.03 eV/A
        from belt4's cell is 0.18 kcal/mol per monomer above the same phase at 0.004, and the 0.02 cell of
        3.14 (`/opt/dft/rows14an`) is 0.04 per monomer (1.15 per cell) above it. `/opt/dft/registry14` (the
        all-trans cells at other shears, 0.03) decides nothing.
        THE TWO CHEAPEST SHEETS IN THE POTENTIAL AT 0.004 eV/A (09:17 UTC; `docs/predictions/
        sheets14_potential_sheared.json`, cells in `mlip/sheets14_sheared_in_all_s1.extxyz`, test
        `test_three_evaluators_put_the_leaders_cheapest_sheet_three_kcal_above_the_film_cells`): FILM CELL
        site 12 `@0` -3.60 (belt4 -2.27, chain1-cfe -1.65); LEADER site 10 `@0` +0.05 (+0.65, +1.66). The gap
        is 2.9 to 3.7 in all three. How: `row_cells_frames.py kink --scan docs/predictions/
        kink_site_scan_belt4.json --key "ROW@0" --sites N --sheared docs/predictions/
        kink_sheet_shear_belt4_sites.json --kmesh "4 2 1" --out /opt/dft/sheets14s`; `python -I
        examples/mlip_train.py relax --model docs/predictions/mlip/models/all_s1.model --name all_s1 --cells
        START --frame NAME_rN --out /opt/dft/sheets14s --threads 2 --fmax 0.004 --steps 300` in legs (each
        leg starts from the last one's cell; half an hour to an hour a leg); `row_cells_frames.py
        sheets-potential --pairs "ROW@SITE=ALL_TRANS.start.xyz,SHEET.start.xyz" ... --out
        docs/predictions/sheets14_potential_sheared.json` (give every pair again: the file is written whole).
        All-trans references: film cell `/opt/dft/rows14/row_vdf9-trfe4-cfe1-14_TTTT_chain1_in_all_s1`
        (0.005), leader `/opt/dft/sheets14s/tt_leader_from_rows14an_r1` (0.004).
        THE POTENTIAL AT 0.004 eV/A ON ALL FOUR CELLS (11:56 UTC 2026-11-20; `docs/predictions/
        sheets14_potential_sheared.json`, `rows14_potential_tight.json`, twelve cells in `mlip/
        sheets14_sheared_in_all_s1.extxyz`, `four_axes_potential_tight_margins_shear.json`; docs 3.20; tests
        in `tests/test_four_axes_at_margins.py`). CHEAPEST SHEET (belt4's cheapest site, one a cell), kcal/mol
        per kink pair: FILM CELL -3.60 (belt4 -2.27), LEADER +0.05 (+0.65), NITRILE ALONE +0.36 (-2.08),
        MATRIX +3.13 (+0.94). So belt4's "nitrile cell with the film" and "leader with the matrix" do not hold
        in the potential; "film cell first, leader three behind" holds in all three evaluators.
        MARGINS (lowest helix placement only): MATRIX +0.016 (3.14: -0.024), FILM CELL +0.129 (+0.078),
        NITRILE ALONE -0.109 (-0.118), LEADER -0.070 (-0.055). At +0.129 the film cell's row reads tan delta
        0.061, 0.071%, 6.2 pN/nm2, on the measured film's channel.
        IN THE POTENTIAL THE FILM CELL'S LOWEST STATE IS THE KINKED ONE: per chain of fourteen against
        all-trans, helix -1.81 and the sheet at the chlorine -3.60 (matrix -0.23, +3.13; leader +0.98, +0.05;
        nitrile +1.53, +0.36). Both presets have the helix lower (belt4 -4.4 against -2.3; chain1-cfe -5.0
        against -4.5). If the potential is right the margin as defined does not name the film cell's lowest
        state. THE FIRST QUESTION FOR DFT AT THIS SIZE.
        How the tight cells were made: copies in `/opt/dft/sheets14s` (`helix_ROW_CONF`, `tt_ROW_...`,
        `sheet_..._sheared`), `mlip_train.py relax ... --fmax 0.004 --steps 300` in legs (`NAME_rN`); margins by
        `row_cells_frames.py potential --frames /opt/dft/rows14tight --rows ...` from cells copied there
        under the names `row_ROW_CONF_TAG_in_all_s1.start.xyz`; sheets by `sheets-potential --pairs` (every
        pair given again).
        THE LEADER'S `@0~` SITE 18 IN THE POTENTIAL IS NOT AT REST AND IS ALREADY BELOW +0.05 (14:03 UTC):
        from belt4's upright sheet (`sheet_vdf9-trfe3-an1-rvdf1-14_belt4_at18m_sheared`, legs `_rN`), per
        kink pair after legs 1 to 5 (1,300 steps, none converged, largest force 0.01 to 0.02 eV/A, kink
        kept): +2.16, +1.01, +0.38, +0.07, -0.28. The cell is sliding: the chain axis tilts steadily (c along
        a and b: -0.12 -0.10, -0.20 -0.20, -0.32 -0.28, -0.37 -0.31, -0.47 -0.39 of a and b), which is a shear
        along the chains of -0.4 and -0.7 monomers so far. SO +0.05 IS NOT THE LEADER'S FLOOR IN THE
        POTENTIAL, and the 3.65 gap to the film cell is an upper bound until this stops.
        16:03 UTC: THE CONTROL IS IN AND THE REFERENCE DOES NOT MOVE. The leader's all-trans cell given the
        sheet's tilt (`/opt/dft/sheets14s/tt_leader_tilted_like_at18m`, legs `_r1` to `_r3`; built by the
        scratch script `tt_tilted_like.py FLAT LIKE OUT FRAME`) converged at 0.004 eV/A still tilted (c along
        a, b: -0.42, -0.27) at -0.035 kcal/mol per chain of fourteen against the upright cell. THE SHEET after
        legs 6 and 7: -0.67, -1.00, not converged, kink kept, tilt -0.55 -0.50, shear along the chains -0.45
        -0.93 monomers. Interim paragraphs are in docs 3.20 and the head of CANDIDATE_MATERIALS; the stored
        JSON, the extxyz and the tests still carry only site 10 (+0.05) and have to be extended when the
        sheet has stopped.
        17:36 UTC: THE SLIDE HAS AN END NEAR A SHEAR OF -0.5 AND -1.0 TO -1.1 MONOMERS. Starts pushed along
        the slide (scratch script `slide_on.py PREV NOW OUT FRAME K`: the chain axis of NOW moved K times
        further than it moved from PREV to NOW, repeat length and fractional coordinates kept; from legs 5 and
        7) come BACK: `..._at18m_slide4` from shear -0.61 -1.33 to -0.51 -1.14 in two legs (-1.02, -1.32
        kcal/mol per kink pair), `_slide8` from -0.75 -1.69 to -0.54 -1.23 (+1.89, -0.86), while the plain
        sheet's leg 8 goes on forward to -0.47 -0.97 (-1.24). None converged; all three still falling about
        0.3 a leg, kink kept in all. A leg of 250 steps takes 18 minutes on four cores (DO NOT RUN TWO AT ONCE
        HERE: side by side on two cores each, neither finished a leg in 38 minutes and both wrote nothing).
        RUNNING since 17:40 UTC, job `b0b2du9vx`: `_slide4` alone, three legs of 500 steps (`STEPS=500
        TMO=2600 legs.sh NAME NAME NAME`: each leg from that name's latest stored `_rN`; a name whose log says
        `converged True` is skipped), so the optimizer keeps its Hessian longer.
        If either is cut, go on from its last `_rN` cell. Read a leg with `row_cells_frames.py
        sheets-potential --pairs "KEY=/opt/dft/sheets14s/tt_leader_from_rows14an_r1.start.xyz,LEG.start.xyz"
        --name all_s1 --all-trans docs/predictions/rows14_potential_tight.json --out SCRATCH.json`.
        THEN: at scale it is A POD'S JOB (108 sheets, 16 more helix cells, tight), with DFT on the film
        cell's all-trans, helix and kinked cells first (168 atoms), then the leader's.
        DFT AT FOURTEEN MONOMERS DOES NOT FIT THIS MACHINE: the leader's helix cell (170 atoms, 350 bands, 6
        k-points) asks for 15.7 GB (pw.x's own estimate at 3 processes, one pool) and there are 15; STOPPED
        before it ran out (`/opt/dft/rows14dft` holds the two cells and the anchor check). DO NOT START IT HERE
        AGAIN. It is a rented machine's job (more memory), which needs the PI's go.
        THE OTHER UNITS (TFE, CDFE, ethylene, VF, DFE): UNTESTED UNTIL THEIR CHECK. Rows by unit: AN in 7 of the
        35, CDFE 4, TFE 2, ethylene 2, VF 1, DFE 1 (15 rows in all); CFE 9, CTFE 3; 9 rows VDF, TrFE and reversed
        VDF only. `examples/unit_dft_check.py --row ROW --work /opt/dft/unit4 --model
        docs/predictions/mlip/models/all_s1.model --pw /opt/mm/qe/bin/pw.x --pseudo /opt/pseudo --check
        /opt/dft/chainsplit/check --threads 3 --np 3` runs the whole check for one small row and goes on from
        where it stopped (about an hour a unit). Rows `vdf2-trfe1-{cdfe,tfe,vf,e,dfe,rvdf}1-4` are in
        `candidate_screen.SEQUENCES`. CDFE (docs 3.13): VRVD DFT +0.968 = +2.394 inside -1.426 between; chain1
        +0.975 = +3.091 -2.117 shear-free: THE MARGIN TO 0.007 BY TWO ERRORS OF 0.7 THAT CANCEL, confirmed by the
        lone VDF-CDFE chains (`chain_conformers_cdfe_pbe_d3bj.json`: helix -2.843 and +3.750; chain1 -3.418 and
        +1.467). NOT COUNTED AS PASSED. TFE (docs 3.13): VRVF DFT -0.076, chain1 -0.060 (parts +0.26, -0.29 off):
        PASSES; its row `vdf9-trfe4-tfe1-14` is +0.147 under chain1-cfe. ETHYLENE AND VF: polyfind's side done
        (chain1 -0.232 and +0.479; belt4 +0.500 and +0.452), their cells relaxed in the potential (VF four of
        five); DFE not started; NO DFT ON THE THREE YET (a container restart at about 13:00 UTC stopped it):
        `unit_dft_check.py --row vdf2-trfe1-{e,vf,dfe}1-4 ...` goes on from there. `--no-dft` stops a check after
        the potential's relaxations. THEN: e, vf, dfe. THE MACHINE HAS FOUR CORES:
        DFT ON THREE (`--np 3`) AND ONE OTHER JOB, NO MORE (a fifth process doubled a DFT step's time). A unit
        that fails INSIDE THE CHAIN gets lone-chain conformers (`chain_conformers.CHAINS` has an2, vf2, cdfe2,
        e2, tfe2, dfe2, ctfe2: all-trans and the helix placements that differ; `frames`, `settle`, `collect`),
        then `unit_torsion_fit.py --base pvdf-dft-chain1-cfe --type ITS_BOND`; one that fails BETWEEN CHAINS gets
        `unit_nonbonded_fit.py`; either is accepted on a cell that was not in the fit.
        THE 35 ROWS, MARGINS UNDER `pvdf-dft-chain1-cfe`: DONE FOR 25 ROWS (2026-11-19 12:50 UTC):
        `docs/predictions/margins_chain1-cfe.json`, table in docs 3.13. The 20 rows whose units pass (TrFE,
        reversed VDF, CFE, CTFE at the edge), the TFE row, and the 4 CDFE rows (unit NOT passed: read them with
        that). In the window of 55 V/um on this axis: `vdf8-trfe4-ctfe1-rvdf1-14` +0.061, `vdf6-trfe4` +0.083,
        `vdf9-trfe4-cdfe1-14` +0.089, `vdf9-trfe4-tfe1-14` +0.147. Below zero: `vdf5-rvdf1-6` -0.199, `vdf9-trfe4-rvdf1-14` -0.164,
        `vdf7-trfe2-cdfe1` -0.095, `vdf8-trfe4-cfe1-rvdf1-14` -0.070, `vdf12-rvdf2-14` -0.030. Matrix
        `vdf10-trfe4-14` +0.222, terpolymer `vdf9-trfe4-cfe1-14` +0.359, `vdf13-trfe6-cfe1-20` +0.307. THE 10 ROWS
        NOT RUN: 7 with a nitrile (no preset has it), and those with ethylene, VF, DFE (unchecked). MARGINS
        ARE THE STRAIN AXIS ONLY: tan delta, frequency and blocked stress per area are NOT redone under this
        preset (the kink tools `dilute_kink.py`, `kink_migration.py`, `kink_site_scan.py` take `--preset`; the
        CFCl bond's barrier is not fitted). How it was run: `python examples/rows_under_preset.py --preset
        pvdf-dft-chain1-cfe --redo-failed --minutes 100 --rows ROW,... [--part K/N] --out
        /opt/dft/margins_chain1-cfe.json` (the --out file is the checkpoint, shared under a lock). The
        checkpoint was seeded with the six rows of chain1's run that hold no CFCl carbon (the two presets are the
        same potential for them; marked `carried_from`). chain1's own run (`/opt/dft/margins_chain1.json`, 14
        rows) is kept for the comparison; its chlorine rows are superseded. The 15 rows with an unchecked unit
        wait for their unit.
        ALSO WAITING: one more DFT cell for VRVC: its TGTG' helix placement is 0.015 below G'TGT in the potential
        and 0.217 above it in chain1-cfe (`/opt/dft/cfe4/row_vdf2-trfe1-cfe1-4_TGTGp_chain1_in_all_s1.start.xyz`,
        not labelled). The fourteen-monomer cells in the potential (`/opt/dft/rows14`: both all-trans cells
        relaxed, the helix cells not; 168 atoms, several seconds a step) were stopped for the DFT.
        THE 35 ROWS UNDER chain1 ALONE (SUPERSEDED FOR THE CHLORINE ROWS; how it was run): three side by side (`for k in 1 2 3; do python
        examples/rows_under_preset.py --preset pvdf-dft-chain1 --redo-failed --minutes 105 --part $k/3 --out
        /opt/dft/margins_chain1.json > LOG_$k & done; wait`, one thread each; the --out file is the checkpoint,
        written under a lock so the parts share it; run the same command again after every stop or restart; about
        7 minutes a cell on average, the slow ones the cells that need the trust box; 168 cells; logs
        `/opt/dft/margins_chain1*.log`). 12 of 35 rows had a margin on 2026-11-18 23:50 UTC: `vdf5-rvdf1-6` -0.199
        (upright +0.042, strain +9.0%: the acceptance's numbers again), `vdf9-trfe4-cfe1-14` +0.259 (upright +0.153;
        belt4 +0.317), `vdf7-trfe2-cfe1` +0.400, `vdf7-trfe2-cdfe1` -0.095, `vdf6-trfe4` +0.083.
        Copy the file to `docs/predictions/margins_chain1.json` when all 35 rows have a margin. The kink tools
        (`dilute_kink.py`, `kink_migration.py`, `kink_site_scan.py`, `polish_margins.py`) take `--preset`; none has
        been run under chain1 yet.
        STILL TO DO, IN ORDER: (1) THE 35 ROWS UNDER chain1 ON ALL FOUR AXES (re-polish each row's cells upright and
        with the shear free, pendant tilts; then strain, tan delta, frequency, blocked stress per area, which need
        the slide and kink barriers redone under chain1 and the chain stiffness of the line); E5 under it (the
        registry rows are 0.77-0.88 of DFT, inside the factor 1.5); (2) a joint polish of all 27 variables with
        every row (hh lines and candidate rows included), and the reversed unit's own angles and bond lengths;
        (3) hh4 conformers (24 atoms, frames built in `/opt/dft/chainconf`: hh4_ff, hh4_hh) to tell the CF2-CF2
        bond from the CH2-CH2 bond, which every row so far turns together. SimpleFF (oligomer/RIS side) does not
        know the two charge terms. Text: `docs/CRYSTAL_STAGE.md` 5.24 (12-6 cannot), 5.25 (the wall), 5.26 (the
        charges, the rigid cells, the candidate's line), 5.27 (the conformers, the torsion, stage 13), 5.28 (the two
        runs, cohesion, the acceptance's two outcomes), 5.29 (the inside of a chain, the line, chain1).

**A machine-learned potential as a second opinion: pretrained ones fail; our own, fitted long on a rented card, is the
DFT to 0.01 kcal/mol where it has frames near and off by 0.07 to 0.21 on a packing it was not shown; relax in it, then
one DFT single point, 2026-11-16 (PI: "Go on 1", then "You're using ots solutions? Thought we
would train our own?").** Text: `docs/MLIP_SECOND_OPINION.md` (3.3 has the table, 3.5 what is needed).
THE TEST: PVDF's beta minus alpha, +0.726 kcal/mol per monomer in PBE-D3(BJ), the one margin held from DFT; a second
opinion has to resolve 0.1 (the switching window is 0.12 wide).
PRETRAINED, DONE, NONE PASSES (`examples/mlip_second_opinion.py`, `docs/predictions/mlip_second_opinion.json`,
`tests/test_mlip_second_opinion.py`): MACE-MPA-0, MP-0b3, MP-0a (MIT, with D3(BJ) added) and Egret-1, Egret-1T (MIT)
are off by 1.2 to 3.9, outside the 0.62 to 1.55 that published functionals span, no agreement on the sign. NOT RUN:
models under academic or gated licences (MACE-OFF23, OMAT, MATPES, OMOL, UMA): the PI's decision.
OUR OWN, DONE SO FAR:
(1) the anchor's DFT recipe runs in the cloud session (QE 7.5, the anchor's energy to the last digit, 52 s).
(2) PVDF frames, 76 (ends 4, strains 24, rattles 18, dynamics snapshots 30; 2.5 h of DFT on four cores):
    `docs/predictions/mlip/pvdf_frames_pbe_d3bj.extxyz`. `examples/mlip_frames.py make / label / collect / cutoff`.
(3) two fits on from MACE-MP-0b3, `examples/mlip_train.py prepare / train / evaluate / closure`, results in
    `docs/predictions/mlip_own_pvdf_ft_0b3.json` (60 epochs, 1 h) and `mlip_own_pvdf_ft_0b3_b.json` (30 more with energy
    and stress weighted up, 37 min); `tests/test_mlip_own.py`. beta minus alpha: +0.671 / +0.700 relaxed (first fit),
    +0.617 / +0.566 (second); DFT +0.726. Antipolar minus polar beta +2.7 and +4.2 meV (DFT -1.1). C_11 39 and 31 GPa
    (DFT 22.8); relaxed volumes 0 to 7% off. DFT single points on the cells each fit relaxed
    (`docs/predictions/mlip/pvdf_probe_frames_pbe_d3bj.extxyz`): forces left 0.01 to 0.08 eV/A; gaps between cells off
    by 0.011 to 0.246 kcal/mol; the two fits differ from each other by more than 0.1 on every gap.
    READING: ten times better than pretrained, at the edge of the resolution needed and not inside it; atoms right,
    lattice across the chains not; the fits are under a thousand optimiser steps each (4 s a step on four cores).
    BY-PRODUCT: DFT puts the helix cell with one chain turned half a turn (dipoles adding, delta-PVDF's kind) level
    with alpha (0.4 meV per monomer below alpha's minimum on one geometry, 2.0 above on the other, neither relaxed in
    DFT). The anchor's +0.726 is against alpha; against the lowest helix packing found it is +0.735 at least.
NOT DONE, IN ORDER:
(a) THE LONG FITS ON A RENTED CARD ARE DONE (PI: "Go on runpod", 2026-11-16). docs/MLIP_SECOND_OPINION.md 3.6.
    Pod `wylsqwj1no08y4`, Runpod Secure Cloud, one A40 at $0.59/h, 17:02:34 to 20:44:31 UTC (3.7 h, about $2.18;
    spend now about $16.23 of $46). It removed itself after `stop` (its status page went from 200 to 404 at
    20:44:33). THE REMOVAL IS NOT CONFIRMED BY AN API READ: since a container restart at about 18:50 UTC every
    request to `rest.runpod.io` answers 502 "upstream request failed" through the session's proxy (the key is still
    on the environment; other hosts answer). THE PI LOOKED (2026-11-16, their words: "The pod has run to
    completion."), so it is not running and not billing for the card. Whether it is deleted or only listed as
    exited was asked and not yet answered. Retry
    `curl https://rest.runpod.io/v1/pods/wylsqwj1no08y4` (it must say the pod does not exist); if the pod is still
    listed, DELETE it and read again.
    WHAT RAN: `examples/gpu_jobs/mlip_fit.sh` (queue `queue_20261116_mlip.txt`; its own Python 3.12 environment,
    torch 2.5.1 for CUDA 12.1): four fits of 400 epochs on the 130 PVDF frames (76 + 54), stress weight 100 or 1000,
    two seeds: `docs/predictions/mlip_own_pvdf_s1.json`, `_s2`, `mlip_own_pvdfb_s1.json`, `_s2`.
    RESULT: beta minus alpha +0.721 to +0.732 on the DFT geometries (DFT +0.726); antipolar minus polar beta -0.5 to
    -1.2 meV (DFT -1.1); C_11 21.4 to 23.7 GPa (DFT 22.8), C_12 40% high. DFT single points on the cells each fit
    relaxed: gaps on the trained packings right to 0.014 kcal/mol (one fit 0.083), THE UNSEEN TURNED-CHAIN HELIX
    PACKING OFF BY 0.07 TO 0.21, always too low; DFT on the four fits' turned cells agrees with itself to 1 meV.
    Relaxed strictly (fmax 0.002, `mlip_own_pvdfb_s1_tight.json`, `_s2_tight`): cells within 0.8% of DFT's volume,
    DFT energy within 0.5 meV per monomer of its own minima.
    READING, AND HOW IT IS USED: relax a cell in the potential, then ONE DFT single point for its energy. The
    potential's own energy is trusted only for kinds of cell it has been checked on.
    FILES OUTSIDE THE REPOSITORY (lost with the container): the four models `/opt/dft/gpu13/*.model` (13 MB each;
    saved on the card: load with `mlip_second_opinion.file_calculator`), DFT single points `/opt/dft/probe13`,
    `/opt/dft/probe13t` (all 32 are in `docs/predictions/mlip/pvdf_probe_frames_pbe_d3bj.extxyz`).
    LESSONS FOR THE NEXT POD: one fit at a time on a card (four side by side took four times as long each); 7 s an
    epoch for 112 small frames on an A40; a job must end well inside the pod's own limit; fetch results as they
    appear (`scratchpad/runpod/mirror13.sh` did). The image's Python 3.10 is too old for mace-torch 0.3.17.
    NEXT (docs 3.7): when the units frames are labelled, one fit on everything on a second pod (ASK THE PI FIRST:
    a new spend), then a candidate's own cells relaxed in it and one DFT single point each.
(b) FRAMES FOR THE OTHER UNITS: DONE, 2026-11-17 04:01 UTC. 168 frames in 36 cells (TrFE, CFE, CTFE, the nitrile,
    reversed VDF, TFE in the PVDF cells, all-trans and helix; 24 to 49 atoms), 13.5 h of DFT on four cores:
    `docs/predictions/mlip/units_frames_pbe_d3bj.extxyz`. `mlip_frames.py make --set units`, `label --only units`,
    `aim` (the sampler's volumes were wrong both ways; the other 132 frames were moved by the DFT stress on each
    relaxed cell), `label`. In all 36 cells the labelled frames have the DFT stress on both sides of zero along both
    lateral axes. The cutoff check passes (`docs/predictions/mlip_cutoff_check.json`: a CFE pair and a nitrile pair
    at 100/800 Ry, gaps change by 0.002 kcal/mol). docs 3.4.
(c) PVDF FRAMES NEAR THE MINIMA: DONE, 54 (`docs/predictions/mlip/pvdf_near_frames_pbe_d3bj.extxyz`), in the long fits.
    THE WHOLE SET IS 330 FRAMES: `python -I examples/mlip_train.py prepare --frames` the four files `--data DIR
    --hold-out "VAVF|VRVV"` gives 282 in the fit, 48 out (one composition whole). It was run here and passes.
    THE PI SAID "Go" (2026-11-17) AND SAVED THE KEY AGAIN (the secure form, 19:07 UTC): `rest.runpod.io` STILL
    ANSWERS 502 in 0.3 s from this session, so no pod can be made from it. The pod needs a NEW SESSION (asked of the PI).
    A SECOND ROUTE THAT DOES NOT USE THE KEY OR THAT HOST (PI asked 2026-11-18, "U using the runpod mcp?"; the answer
    was no: every pod so far was made with curl against the REST API): Runpod's hosted MCP server,
    `https://mcp.runpod.io/` (https://agents.runpod.io/; sign-in by OAuth, no API key; it answers from the session).
    It is NOT connected: no Runpod connector on the PI's claude.ai account, no Runpod tools in the session; the plugin
    installed on 2026-11-13 lived in that container only and its sign-in was never done. For a cloud session the
    PI adds it as a connector at https://claude.ai/customize/connectors (URL `https://mcp.runpod.io/`, sign in with
    Runpod) and starts a NEW session: connectors are read when a session starts. Its tools are untested here: check
    that a pod can be made with an image, a start command, environment variables and an HTTP port before relying
    on it; if not, it can still read and delete pods.
    IN A NEW SESSION: check `curl -sS -o /dev/null -w "%{http_code}" https://rest.runpod.io/v1/pods` gives 200; confirm
    stock and price (public GraphQL, as in docs 3.5); state the hourly cost; POST https://rest.runpod.io/v1/pods with
    name polyfind-mlip-fit-all, imageName runpod/pytorch:2.2.0-py3.10-cuda12.1.1-devel-ubuntu22.04, gpuTypeIds
    ["NVIDIA A40", "NVIDIA RTX A6000"], gpuCount 1, computeType GPU, cloudType SECURE, interruptible false,
    containerDiskInGb 40, volumeInGb 0, ports ["8000/http"], env BENCH_TOKEN (32 hex, kept in the scratchpad, never
    shown or committed), BENCH_QUEUE examples/gpu_jobs/queue_20261117_mlip.txt, BENCH_LIMIT_S 18000, BENCH_IDLE_S
    2400, BENCH_LINGER_S 3600, MLIP_HOLD_OUT "VAVF|VRVV", dockerStartCmd ["bash","-c","curl -sL
    https://raw.githubusercontent.com/Goober68/polyfind/<head commit>/examples/runpod_jobs_pod.sh | bash"] (do NOT
    send allowedCudaVersions); read the pod back; after the timing job add ONE line to the queue for the long fit
    (`mlip_fit.sh all EPOCHS 1 frames <the four files> train --stress-weight 1000 --energy-weight 10000`, EPOCHS so
    that it ends well inside the limit) and push; fetch results as they appear (STATUS, jobs/NN_mlip_fit/*.model,
    mlip_own_*.json, probe/, log.txt at https://<pod>-8000.proxy.runpod.net/<token>/); push `stop`; confirm by a read
    that the pod is gone. The model (13 MB) then has to reach the machine with the DFT program.
    THE PLAN AS IT WAS: a second pod, queue `examples/gpu_jobs/queue_20261117_mlip.txt`
    (a timing job first, then ONE fit of everything, stress 1000 / energy 10000; set the pod's `MLIP_HOLD_OUT` to
    `VAVF|VRVV`). AND IT NEEDS `rest.runpod.io` TO ANSWER FROM THE SESSION: still 502 after three container restarts
    (only that host, the one with the saved key); the PI saving the key again through the secure form
    (`request_credential`) is the first thing to try, a new session the second.
    THE CONTAINER RESTARTED AT 18:50, 00:04 AND 01:47 UTC: files kept, the running job killed each time.
(d) THE TEST THAT MATTERS, FIRST PASS DONE, NOT SETTLED (2026-11-17, docs 3.8, `docs/predictions/
    mlip_candidate_vdf5-rvdf1-6.json`, frames `docs/predictions/mlip/candidate_frames_pbe_d3bj.extxyz`). PVDF with one
    reversed unit in six, 72 atoms, polyfind's margin +0.064. DFT single points (15 to 18 min each here):
    on polyfind's cells as they are, all-trans minus helix +0.164 kcal/mol with 1.6 to 2.3 eV/A left on the atoms;
    on the same cells relaxed in the PVDF-only potential `pvdfb_s2` (`mlip_train.py relax`), -0.276 with 0.6 to 0.9
    eV/A left. Relaxing lowered DFT's energy by 43 (all-trans) and 24 (helix) meV per monomer.
    READINGS: (1) a DFT single point on polyfind's cell does not test polyfind's margin at the 0.1 level (this is
    how E7's first-step answers from Dipole must be read); (2) the PVDF-only potential is outside its frames on
    reversed units (0.57 kcal/mol rms between frames of the same atoms on the 56 H/C/F frames of the units set), so
    -0.276 is not the DFT margin; (3) THE SIGN OF THIS CANDIDATE'S DFT MARGIN IS OPEN.
    RELAXING ON IN DFT ITSELF (`mlip_frames.py relax --frames /opt/dft/cand_dft --names rvdf6_TT,rvdf6_TG --fmax 0.05
    --minutes 100`, runs that end by themselves, about 20 min a step; frames in `docs/predictions/mlip/
    candidate_relaxation_frames_pbe_d3bj.extxyz`; the step table is in the result JSON and docs 3.8). AFTER 5 TO 6
    STEPS EACH (2026-11-17 19:10 UTC, run 3 started then): forces 0.06 to 0.11 eV/A, the gap -0.276, -0.281, -0.271,
    -0.269, -0.268, -0.267. SO DFT PUTS ALL-TRANS 0.27 BELOW THE HELIX WHERE POLYFIND PUTS IT 0.06 ABOVE: 0.33 apart,
    6.0 against 4.1 kcal/mol per reversed unit. NOT FINAL (target 0.05 eV/A) AND NOT A SEARCH: the helix cell's other
    arrangements (axial registry) were not looked for. If it holds, every cell placed by a reversed unit sits about
    0.13 per monomer lower at one in fourteen. FINISHED 20:50 UTC: BOTH RELAXED IN DFT TO 0.05 eV/A (7 and 10 steps;
    17 frames, 5.3 h), THE MARGIN -0.261 (polyfind +0.064; 5.9 against 4.1 per reversed unit).
    OTHER HELIX ARRANGEMENTS, DONE 2026-11-18 00:47 UTC: the DFT-relaxed helix cell with its neighbours along b or a
    moved by whole monomers along the chain (seven cells; each chain made whole before the lattice vector is
    changed, or the chains break) and polyfind's other helix conformer GTG'T, each relaxed in `pvdfb_s2` and given
    one DFT single point (`docs/predictions/mlip/candidate_registry_frames_pbe_d3bj.extxyz`). Against the helix cell at
    the same stage: along b -0.9, -1.2, +0.3, +0.2, -0.4 meV per monomer; along a +33, +32; GTG'T +8.5. NONE LOWER BY
    MORE THAN 1.2 meV: THE MARGIN IS -0.26 (-0.23 if the lowest holds), POLYFIND +0.064. Not tried: both lattice
    vectors at once, other setting angles, other all-trans arrangements (those can only widen the gap).
    WHERE THE 0.33 COMES FROM, DONE 2026-11-18 01:08 UTC (docs 3.8, last part; `all_trans_cell_held_upright` in the
    JSON): the all-trans cell relaxed in `pvdfb_s2` with the shear along the chains forbidden (`mlip_train.py relax
    --mask "1 1 1 0 0 1"`, as polyfind's cell is: every reversed unit level with its images' along a and b), one DFT
    single point: ALL-TRANS MINUS HELIX +0.070 (POLYFIND +0.064), and 15.0 meV per monomer (0.346 kcal/mol, 2.1 per
    reversed unit) above the sheared cell at the same stage. SO POLYFIND'S ENERGIES ARE RIGHT AND ITS CELL LACKS A
    FREEDOM: the lattice vectors a and b cannot take a component along the chain, so a pin always sits in a column
    of pins 5 A apart. In a random copolymer a pin's neighbours are ordinary units, so the sheared cell is argued
    (not computed) to be the nearer to the material. The column is in every polyfind cell with one pin per chain
    period; its cost per pin does not thin with the period (about 0.15 per monomer at one in fourteen if 2.1 per
    unit carries over); for CFE, CTFE, the nitrile and TrFE it is not measured and need not have the same sign.
    PROPOSED TO THE PI, NOT STARTED: give polyfind's cell the two shear freedoms and polish the candidates again
    (the number to meet: 0.35 per monomer for this all-trans cell, nothing for its helix); and the same DFT pair
    (upright against sheared) for one CFE and one nitrile cell once a potential can relax them.
    WHAT THIS MEANS FOR POLYFIND (written before the upright test; the PI has it, 2026-11-18): a reversed unit lowers the margin by
    5.9 kcal/mol in DFT against 4.1 in polyfind at one in six; if it scales with their number, cells placed by a
    reversed unit sit about 0.13 per monomer lower at one in fourteen than polyfind has them. Possible causes to
    look at: polyfind's rigid PVDF bonds and angles at the CF2-CF2 and CH2-CH2 junctions (DFT finds 2.3 eV/A there),
    and its cell not letting neighbouring chains' defects sit off level (the all-trans cell sheared to do so).
    THE FASTER ROUTE, for this and the larger candidates: the fit on all 330 frames (54 hold a reversed unit), then
    `python -I examples/mlip_train.py relax --model NEW.model --name NAME --cells /opt/dft/cand/cand_vdf5-rvdf1-6_*.dft.xyz
    --out DIR`, `mlip_frames.py label --frames DIR`, and the gap from the two DFT energies; the DFT force left says
    whether to believe it. If /opt/dft/cand is gone: the start files are made from
    `docs/predictions/polished_vdf5-rvdf1-6_{TTTT,TGTGp}_tilted.xyz` with n_monomers 12, kmesh "4 2 2".
RUNNING THINGS HERE: one heavy job at a time (a fit, a DFT run and a sampler run each want all four cores; two
together slow each other by more than half). A background job is stopped at two hours and is not to be relaunched
as it was: use runs that end by themselves (`--minutes`), or foreground runs under ten minutes.
TO RESTORE THE DFT SIDE IN A NEW CONTAINER (nothing of it is in the repository but the scripts and the frames):
    pip install torch --index-url https://download.pytorch.org/whl/cpu; pip install mace-torch ase torch-dftd
    micromamba (https://micro.mamba.pm/api/micromamba/linux-64/latest) create -p /opt/mm/qe -c conda-forge "qe=7.5"
        (run it in the foreground; Ubuntu's own quantum-espresso 6.7 package crashes on reading a pseudopotential)
    SSSP_1.3.0_PBE_precision.tar.gz from https://archive.materialscloud.org/records/rcyfm-68h65 -> /opt/pseudo
        (H_ONCV_PBE-1.0.oncvpsp.upf, C.pbe-n-kjpaw_psl.1.0.0.UPF, F.oncvpsp.upf, N.oncvpsp.upf,
        Cl.pbe-n-rrkjus_psl.1.0.0.UPF; checksums of the first three in deliverables/polar_margin_dft/raw_stages.json)
    python -I examples/mlip_train.py prepare --frames docs/predictions/mlip/pvdf_frames_pbe_d3bj.extxyz --data /opt/dft/train
    python -I examples/mlip_train.py train --data /opt/dft/train --name pvdf_ft_0b3      (the first fit again: 1 h)
    The model files are not in the repository (13 MB each, and a model file is a pickled network).
Egret's model files: `add_repo rowansci/egret-public` (read), clone to /home/user/rowansci/egret-public; a model file
is a pickled network, run the Python that loads one with `-I`.

**THE PI IS THE BENCH (2026-11-15: "Im the bench, so go").** Bench asks go to the PI directly, in the conversation,
and are filed in the exchange for the record. Filed: `sarco/materials/ask-polyfind-bench-film-stiffness-v1.md`
(four numbers on the baseline 25 um P(VDF-TrFE-CFE) film, in order: (1) the storage modulus, from the DMA run E4 already
lists, or the small-strain tensile modulus; (2) crystallinity and draw; (3) the blocked force itself against a stiff
stop at the operating point, with the free strain; (4) the regio-defect content, 19F NMR or the supplier's figure).
WHEN THEY COME: (1) into `pin_design.py --young`; (1)+(2) into a two-phase estimate between the film's modulus and the
crystals' (`crystal_stiffness_belt4.json`); (3) calibrates the blocked-force column on the baseline; (4) into the film's
cell (a cell at the measured defect content) and the anchor of the relaxation frequencies. The bench's instruments and
limits are in `sarco/materials/bench-expectations-v1.md` and `sarco/materials/electrostriction/BENCH_E4.md`.

**The crystals' stiffness, and the blocked force through the crystal, 2026-11-15 (PI: "Go").** Text:
`docs/CANDIDATE_MATERIALS.md` 5.6 ("The fourth part: blocked force, per area": rewritten with both readings). STATE: DONE,
CPU ONLY (spend unchanged: about $14.05 of the $46).
WHAT WAS RUN: `examples/crystal_stiffness.py --rows ROW --cells TTTT,helix` (new): the polished two-chain cell strained
along xx, yy, zz and the shear xy, every internal coordinate relaxed at the fixed cell on the closure surface, the
elastic constants by central differences over 33 cells; ten minutes per fourteen-monomer cell. Thirteen cells in
`docs/predictions/crystal_stiffness_belt4.json`. CHECK: beta-PVDF C_11 24.4, C_22 25.3, C_12 5.2 GPa (DFT 22.8, 24.5, 5.9).
C_33 IS NOT USABLE (rigid bonds: 674 GPa for beta-PVDF, DFT 303).
YOUNG'S MODULUS ACROSS THE CHAINS (mean of the two lateral axes, GPa; all-trans / helix):
    beta-PVDF 23.6 | 71/29 matrix 14.1 / 11.1 | 64/29/7 (the film's composition) 13.7 / 11.5 | nitrile alone 10.6 / 10.2
    PVDF with reversed units 10.1 / 12.1 | nitrile with TFE 8.6 / 9.3 | the leader (nitrile + reversed VDF, 21% TrFE) 8.5 / 7.2
BLOCKED FORCE PER AREA, TWO READINGS (pN/nm2, heat budget, 34 kHz; `pin_design.py`, `--stiffness`):
                                             film's assumed modulus     through the crystal across the chains (upper bound)
    the leader, by the move / by the sum      0.28 / 0.34                7.4 / 8.8
    PVDF with reversed units (sum)            0.16                       7.0
    the film's channel                        0.13                       6.5
    the nitrile alone, move / sum             0.18 / 0.19                6.1 / 6.5
    the nitrile with TFE (one site)           0.19                       5.7
READINGS: (1) every pin softens the crystal across the chains and the leader's is the softest; (2) on the crystal's
reading the ranking flattens: the leader stays first by 1.15 to 1.35x over the film's channel (2.1 to 2.5x on strain
alone), the nitrile alone is level with the film, PVDF with reversed units is second; (3) the two readings are 25 to 50x
apart and a real film is between them: the load path (crystallinity, orientation) decides, which is the bench's modulus.
ASSUMPTIONS of the crystal reading: the two states in series by the polar fraction; the switch's areal strain split
evenly between the lateral axes; rigid pendants; the film's channel on its composition's cell.
Test: `tests/test_kink_step.py::test_the_crystals_stiffness_and_the_blocked_force_through_it`.
NEXT, proposed to the PI: (a) ask the bench for the film's modulus at 34 to 66 kHz and its crystallinity (the number
that places a film between the two readings); (b) the CTFE pair on more TrFE; (c) the stiffness of the CTFE cells.

**The fourth part of the objective: blocked force, per area, 2026-11-15 (PI: "Do we have a way to estimate blocked
force? Thats the 4th part of the equation"; then "Actuator has no fixed xsection, so establish an abstract, say force =
N/nm2").** Both are in `CLAUDE.md`. Text: `docs/CANDIDATE_MATERIALS.md` 5.6 ("The fourth part: blocked force, per area")
and a sentence in section 1. STATE: THE COLUMN IS IN; IT IS NOT YET AN AXIS OF ITS OWN.
WHAT IS IN: `pin_design.py` carries `blocked_stress_pa` and `work_j_per_m3` at the heat budget for every cell (and
`--young`, default 0.2 GPa) and its four tables print "blocked force per area, pN/nm2" (1 pN/nm2 = 1 MPa = 1 N/mm2;
plain N/nm2 would be 1e-13-sized numbers). It is `polyfind.loss` `operating_point`'s `blocking_stress_pa` = Young's
modulus times the AC strain; work per cycle against a matched load is a quarter of their product.
    the nitrile with reversed VDF on 21 mol% TrFE: 0.28 (one chain's move) to 0.34 (sites' sum) pN/nm2, 100 to 141 J/m3 per cycle
    the nitrile alone 0.18 to 0.19 | PVDF with reversed units 0.16 | the film's channel 0.13 | CTFE with reversed VDF 0.02
WHAT IT IS NOT: one assumed modulus for every cell, so it ranks as strain ranks. To make it separate: (a) the crystal's
stiffness per candidate, `polyfind.mechanics.elastic_constants` (and its `Actuator` blocking stress), run so far on
beta-PVDF only (E1); (b) the film's modulus, which is the amorphous fraction's and the morphology's: the bench's number
at 34 to 66 kHz or a bound from an assumed crystallinity. The switching is not the limit: a stress moves the margin by
0.0009 kcal/mol per pN/nm2, half the window at 67 pN/nm2.
Test: `tests/test_kink_step.py::test_the_blocked_force_per_area_is_the_films_modulus_times_the_strain`.
NEXT, proposed to the PI: the crystal stiffness of the four leading cells (CPU), then the CTFE pair on more TrFE.

**Every helix conformer polished: five margins moved and the head of the ranking with them, 2026-11-15 (PI: "Go").**
Text: `docs/CANDIDATE_MATERIALS.md` 5.6 (corrected in place, each correction dated; the bullet "Every helix conformer is
polished since 2026-11-15"), 5.4, 5.5, 8.5, header; `docs/E1_E2_PREDICTIONS.md` ("E7 part d, corrected"). STATE: DONE, CPU
ONLY (spend unchanged: about $14.05 of the $46).
THIS BLOCK SUPERSEDES, IN EVERY BLOCK BELOW: that CTFE with a reversed VDF unit is on the switching window, leads, or
carries 0.16 to 0.19% (its margin is +0.245, not +0.076: it is fast and does not switch); the reversed unit's pin of
-4.8 (it is -3.1) and its margin -0.289 (-0.162); that the nitrile pair's pins add; TFE's -1.5 and ethylene's -1.7
(-1.2, -1.4).
WHAT WAS RUN: `polish_margins.py --rows ROW --tg "TGTG',TG'TG,GTG'T,G'TGT"` for the 28 rows with fewer than four helix
conformers (about 75 minutes on four cores). 24 completed; in the four CDFE rows some conformers cannot be built (a
chlorine within bonding distance of a fluorine: `offset_sites` raises), tried one by one.
MOVED (before -> now, kcal/mol per monomer; the new lowest helix):
    vdf8-trfe4-ctfe1-rvdf1-14  +0.076 -> +0.245  GTG'T (TGTG' +0.169 above it)   THE SCREEN'S LEADER: OFF THE WINDOW
    vdf9-trfe4-rvdf1-14        -0.289 -> -0.162  GTG'T                            the reversed unit's pin -4.84 -> -3.06
    vdf6-trfe3-cfe1            +0.338 -> +0.454  G'TGT                            a ten-monomer row, not in the ranking
    vdf9-trfe4-tfe1-14         -0.050 -> -0.026  G'TGT        vdf9-trfe4-e1-14   -0.066 -> -0.045  G'TGT
STOOD: the 64/29/7 cell +0.317, the matrix +0.056, AN +0.197, CTFE +0.300, AN with TFE +0.186, AN with ethylene -0.119,
PVDF with reversed units +0.162, the six-mer +0.064, 70/20/10 AN +0.565, 80/20 +0.229, and the three rows of today
that had all four from the start (+0.037, -0.148, -0.044).
THE RANKING NOW (three axes; the film's channel 364 kHz, tan delta 0.060, 0.066% at the heat budget):
    1. the nitrile with a reversed VDF unit on 21 mol% TrFE: 2.6 MHz (one chain's move) to 6.9 MHz (sites' sum),
       tan delta 0.011 to 0.004, 0.142 to 0.168%
    2. the nitrile with TFE (one site: 1.6 MHz, 0.097%) and the nitrile alone (1.5 to 1.8 MHz, 0.003 to 0.004, 0.091 to 0.097%)
    3. PVDF with reversed units (0.9 to 1.1 MHz, 0.076 to 0.081%)        4. the film's channel
    out: both CTFE cells (5 to 52 MHz, 0.005 and 0.011%: off the window), the model's CFE cell, every sparse pin.
THE REVERSED UNIT'S PULL GOES WITH ITS NEIGHBOUR: per unit in fourteen 0.36 beside CFE, 0.35 beside AN, 0.22 alone in
the matrix, 0.06 beside CTFE; 4.1 per unit in PVDF. Reading 4 of the film (5.3) rests on the CFE figure and stands.
REGISTERED NUMBERS: `e7_polyfind_belt4.json` is KEPT AS SENT (never write `--set e7` again: it would take the new
ground helix); part d's correction is `e7_addendum_polyfind_belt4.json` (`--set e7b`) with
`polished_vdf8-trfe4-ctfe1-rvdf1-14_GTGpT[.xyz, _tilted.xyz]`. Sent to Dipole: sarco main `a0d6ede3` (the request's README and ask carry a dated correction; the two GTG'T files and
the addendum set are in `materials/electrostriction/requests/unit_screen_cells_e7_v1/`; two rows of the regio-defect
note's table corrected).
Tests: `tests/test_kink_step.py::test_every_helix_conformer_is_polished_and_five_margins_moved`; six earlier tests
re-pointed at the corrected numbers.
TRAPS: (1) NEVER RANK A CELL ON FEWER THAN FOUR HELIX CONFORMERS. `polish_margins.py` without `--tg` takes those within 0.3
of the SEARCH's lowest, and the polish moves a helix by 0.4 to 0.9. (2) THE CONTAINER WAS REPLACED at about 05:27 UTC
on 2026-11-15: the working copy came back empty and the scratchpad was lost; nothing unpushed was lost. To restore:
`git fetch --depth 50 origin claude/polymeric-stable-arrangements-uh06b1`, check it out, `pip install -e ".[dev]"`,
clone `Goober68/sarco` beside it. Push before leaving anything running for long.
NEXT, proposed to the PI: (a) bring the CTFE pair back to the window with TrFE (one more TrFE in fourteen, 36 mol%: the
fastest lattice, if the nitrile pair's lever holds there), all four conformers from the start; (b) the nitrile pair
between 21 and 29 mol% TrFE in a twenty-eight-monomer repeat; (c) the film's cell with the measured defect content when
the number comes; (d) Dipole's E7 a.

**One chain's move in the nitrile cell with a reversed unit on 21 mol% TrFE; the note to Dipole, 2026-11-15 (PI: "Go
both").** Text: `docs/CANDIDATE_MATERIALS.md` 5.6 (the block "The reversed unit beside the nitrile, and beside CFE": the
bullet "By one chain's move..."; the moves table) and 8.3 (the cell from two start sites); `docs/E1_E2_PREDICTIONS.md` (E7:
"What E7 a also decides"); job logs and README in `docs/predictions/gpu_jobs_a100_20261115_t/`. STATE: DONE; THE POD IS
DELETED AND CONFIRMED BY A READ (GET 404, 0 pods; 25 minutes, about $0.91; about $14.05 of the $46 used in all). The
fallback check-in is deleted.
(1) THE MOVES, sixteen chains, per-chain metric (queue `examples/gpu_jobs/queue_20261115_t.txt`; rows
`vdf9-trfe3-an1-rvdf1-14@18 4x2 ...` and `@0 4x2 ...` of `kink_migration_belt4.json`; states in `kink_states_belt4.json`),
kcal/mol against the sheet at the start site:
    from the VDF run (offset 18, sheet +1.86): I +8.46, B +6.79, the step 10.35 / 3.56 back (settled), N +7.03,
        sheet to N 10.22 / 3.19 back, L +6.12, NO PATH FROM N TO L (band not converged)
    from the reversed unit (offset 0, sheet +2.60): I NOT A STATE, B +12.85, N +7.98, sheet to N 10.66 / 2.68 back,
        L +4.48, N to L 2.13 above N (2.5 to L per return)
(2) THREE AXES ON BOTH RATE MODELS for this cell (`pin_design.py`; the film's channel: 364 kHz, 0.060, 0.066%):
    by one chain's move (10.22 against the anchor's 11.40: 7.2x): 2.6 MHz, tan delta 0.011 at 34 kHz, 0.142% at the heat budget
    by the sum over kink sites (18.9x):                              6.9 MHz, 0.004, 0.168%
    The models differ by 2.6x here (the nitrile alone 1.2x, the CTFE cell with a reversed unit 10x). SECOND PLACE STANDS:
    by the move the order is CTFE with reversed VDF 5.1 MHz / 0.160%, this cell 2.6 MHz / 0.142%, the nitrile 1.5 MHz / 0.091%.
(3) THE NOTE TO DIPOLE: `sarco/materials/ask-polyfind-dipole-film-regio-defects-v1.md` (sarco main `b9a1aee7`): E7 a now
also decides reading 4 of the film; asked, with no computation: the regio-defect content of the bench's 63/30/7 film.
Dipole is still on E2's fourth geometry (SCF 10 on 2026-11-14); nothing of E5, E7 has started.
CODE: `examples/gpu_jobs/band_row.sh`, `settle.sh`, `between.sh` take `ROW@OFFSET`; `band_row.sh ... ends` runs the end
states only; `pin_design.py` and `candidate_frequency_map.py --cheapest` take, for a row stored from several start sites,
the move from the cheapest sheet (the others are listed beside it); `kink_summary.py` labels a row's start site.
Test: `tests/test_kink_step.py::test_the_nitrile_cell_with_a_reversed_unit_by_one_chains_move`.
TRAPS: (1) A start site needs its polished two-chain kinked cell first: `kink_migration.py --row ROW --kink-at K --nx 1
--ny 1 --mover all --reference-only` (CPU, about ten minutes), committed before the pod starts. (2) Queue `band_row.sh
... ends` first for a new lattice or site and read I: where I is the sheet (+0.000) do not queue the step's bands. The
wide kink is not a state at a reversed unit's site (two cells now). (3) A sixteen-chain dry run on CPU takes minutes per
state: check a script change by its first lines and stop it.
NEXT, proposed to the PI: (a) the film's cell with the measured defect content, once the number comes, and the anchor
re-read on it; (b) the nitrile pair's TrFE content between 21 and 29 mol% in a larger repeat (a twenty-eight-monomer
cell) to place the margin at the window's middle; (c) the remaining helix conformers of the rows polished earlier;
(d) the CTFE cell's move from its cheapest site (the scripts can now start there).

**The reversed unit beside the nitrile and beside CFE, 2026-11-15 (PI: "Go").** Text: `docs/CANDIDATE_MATERIALS.md` 5.6
("The reversed unit beside the nitrile, and beside CFE"), 5.3 (reading 4 of why the film switches), 8.4 (profiles).
STATE: DONE, CPU ONLY (spend unchanged: about $13.15 of the $46).
THREE NEW ROWS (`candidate_screen.py` SEQUENCES; screen file `candidate_screen_belt4_14mer_f.json`; ALL FOUR helix
conformers polished for each, so these margins are not lower bounds; kink sites scanned for both hands):
    vdf8-trfe4-an1-rvdf1-14  IVRVVRVAVVRVVR  margin -0.148 (pins add: -0.149)  sites' sum 13.2x  4.8 MHz  ALL-TRANS SIDE OF ZERO
    vdf9-trfe3-an1-rvdf1-14  IVRVVRVAVVVVVR  margin +0.037 (21 mol% TrFE)      sites' sum 18.9x  6.9 MHz  tan delta 0.004 at 34 kHz
                                             0.130% at 1.8 V/um, 0.168% at the heat budget (34 kHz); dipole 0.449 e.A
    vdf8-trfe4-cfe1-rvdf1-14 IVRVVRVCVVRVVR  margin -0.044 (CFE alone +0.317)  sites' sum 2.8x   1.0 MHz  ALL-TRANS SIDE OF ZERO
READINGS: (1) THE NITRILE WITH A REVERSED UNIT ON 21 MOL% TrFE IS SECOND ON THE THREE AXES, WITHOUT CHLORINE (the film's
channel: 364 kHz, 0.060, 0.066%; CTFE with reversed VDF: 52 MHz by the sum, 5.1 by one chain's move, 0.16 to 0.19%).
Its rate is the sites' sum only: one chain's move is NOT computed for it. (2) Two levers place a margin: one TrFE in
fourteen is 0.18 (0.026 per mol%), one reversed unit in fourteen 0.35 (0.05 per mol%; 0.22 beside CTFE). The window is
0.12 wide. (3) A cell at or under zero gets the degenerate channel from the loss mapping (section 3): a convention, not
a prediction. The nitrile with ethylene (-0.119) was ranked on it and is now left out of the ranking. (4) READING 4 OF
THE FILM (5.3): the model's 64/29/7 cell is head-to-tail only; ordinary PVDF carries 3.5 to 6.0 mol% reversed units
(Bell Labs 1987, found by web search 2026-11-15; no figure found for the terpolymer); at 0.05 per mol% that puts the
cell at +0.14 to +0.01, the window's middle at 4.9 mol%. The same paper: PVDF goes alpha to beta above 11.4 mol%; the
model's PVDF crosses zero at 18 mol% (direction right, two thirds of the strength). (5) THE ANCHOR: every frequency is
364 kHz times a ratio to the defect-free 64/29/7 cell; if the film is the cell with defects, frequencies fall by up to
2.8x (6.9 MHz to 2.4); ratios between candidates do not change. (6) Beside the nitrile and CFE the reversed unit gives
2.7 to 2.8x on the sum over sites; beside CTFE 7.9x.
CODE: `pin_design.py` OTHER_CELLS (a designed cell on another matrix), COMBINATIONS gained the two 29 mol% cells.
Test: `tests/test_kink_step.py::test_the_reversed_unit_beside_the_nitrile_and_beside_cfe`.
TRAPS: (1) `polish_margins.py` polishes only the helix conformers within 0.3 of the SEARCH's lowest; the polish moves a
helix by 0.6 to 0.9, so pass `--tg "TGTG',TG'TG,GTG'T,G'TGT"` to make a margin a number. It did not change these three,
but the rows polished before 2026-11-15 still carry two conformers at most. (2) The screen list is now EIGHT files
(`_14mer_f`), comma-separated in one `--screen` argument.
REPRODUCE: `candidate_screen.py --sequences ROW --conformations "TTTT,TGTG',TG'TG" --out <file>`; `polish_margins.py
--rows ROW --tg "TGTG',TG'TG,GTG'T,G'TGT"`; `kink_migration.py --row ROW --nx 1 --ny 1 --mover all --reference-only`;
`kink_site_scan.py --row ROW [--mirror-kink]`; `candidate_loss_map.py --screen <eight files>`; `kink_summary.py`;
`candidate_frequency_map.py` (both forms); `pin_design.py`.
NEXT, proposed to the PI: (a) one chain's move in the nitrile cell with a reversed unit on 21 mol% TrFE (a pod, under an
hour: the other rate model for the new second place); (b) tell Dipole that E7 a (the reversed unit in PVDF) now also
decides reading 4 of the film, and ask the bench or the literature for the terpolymer's regio-defect content; (c) the
remaining helix conformers of the rows polished earlier (CPU); (d) the CTFE cell's move from its cheapest site.

**The leaders by one chain's move, and CTFE without the reversed unit, 2026-11-15 (PI: "Go on both").** Text:
`docs/CANDIDATE_MATERIALS.md` 5.6 ("The leaders by the other rate model, and CTFE without the reversed unit"), 8.3 ("The
same moves in the unit screen's leaders"), 8.4 (CTFE alone's profile), 8.5; job logs and README in
`docs/predictions/gpu_jobs_a100_20261115_s/`. STATE: DONE; THE POD IS DELETED AND CONFIRMED BY A READ (GET 404, 0 pods;
42 minutes, about $1.53; about $13.15 of the $46 used in all). The fallback check-in is deleted.
(1) ONE CHAIN'S MOVES, sixteen chains, per-chain metric (queue `examples/gpu_jobs/queue_20261114_s.txt`; rows
`<row>@0 4x2 mover 0 ahead host flexible` of `kink_migration_belt4.json` with `to_N` and `to_L`; states N, L and bands A-N,
N-L in `kink_states_belt4.json`), kcal/mol against the sheet:
    AN cell (vdf9-trfe4-an1-14): I +8.87, B +11.86, step 14.98 (band top; its settle slid), N +7.95, sheet to N 10.56 / 2.61
        back, L +4.45, N to L 2.05 above N (2.6 to L per return), L back 6.11
    CTFE with reversed VDF (vdf8-trfe4-ctfe1-rvdf1-14): I NOT A STATE (relaxes to the sheet), B +11.03, NO PATH TO B, N +7.25,
        sheet to N 9.83 / 2.58 back, L +3.40, N to L 2.94 above N (0.54 per return: the saddle to L, 10.19, is the higher), L back 6.79
    for comparison 64/29/7: to N 11.40 / 2.58; 71/29: 11.89 / 2.95
(2) THREE AXES ON BOTH RATE MODELS (`pin_design.py`, tables "sites" and "moves"; frequency, tan delta at 34 kHz, AC strain
at the heat budget at 34 kHz; the film's calibrated channel is 364 kHz, 0.060, 0.066%):
    AN, one per 14:            by the move 1.5 MHz, 0.004, 0.091% | by the sites' sum 1.8 MHz, 0.003, 0.097%   (AGREE)
    CTFE with reversed VDF:    by the move 5.1 MHz, 0.006, 0.160% | by the sites' sum 52 MHz, 0.001, 0.189%   (A DECADE APART)
    71/29 matrix:              by the move 158 kHz, 0.198, 0.054% | by the sites' sum 33 kHz, 0.889, 0.026%
    The move's barrier spans 2.1 across the four lattices where the sheet's energy spans 6.1: that is the two models'
    whole disagreement. Carry the lower figure for the CTFE cell.
(3) CTFE ALONE, one per 14 (`vdf9-trfe4-ctfe1-14@0` and `@0~` in `kink_site_scan_belt4.json`): sites 1.8 to 5.6, wall VK 14.2
for one hand; 18x the anchor's sum (6.6 MHz); margin +0.300, OFF THE WINDOW: 0.1% polar, 0.005% at the heat budget. The
reversed unit brings the margin to +0.076 (the strain) and lowers twelve of thirteen sites by 1.1 to 3.1 (the mirrored hand:
0.7 to 1.7 at ten of twelve), +1.4 at the site before itself; the sum rises 7.9x. THE TWO UNITS' CUTS ADD (0.84 computed, 0.87
summed); THEIR PINS DO NOT (+0.076 computed, -0.046 summed).
CODE: `kink_summary.py` and `candidate_frequency_map.py --cheapest` take a row with a way to N and no step to B;
`pin_design.py --steps` (table "moves"). Test: `tests/test_kink_step.py::test_the_leaders_by_one_chains_cheapest_move_and_ctfe_alone`.
TRAPS: (1) A queue line cannot set the kink offset: both rows start from the sheet at offset 0, which in the CTFE cell with
a reversed unit is the site at the reversed unit (1.00), not the cheapest (0.35). (2) Where the wide kink I is not a
state, `band_row.sh` still runs the "I-B" band, as a straight path from the sheet to B across cis: 900 iterations and
17 minutes of the pod for nothing. Settle I first and look at its energy before queueing the bands of a new lattice.
(3) The AN row's settle of the I-B top slid (to +7.85): stored as `search_left_the_path`; the band's top is the number.
(4) The CTFE row's failed band is stored as `step_to_B_not_found`, without `barrier_forward`.
REPRODUCE: `kink_site_scan.py --row vdf9-trfe4-ctfe1-14 [--mirror-kink]`; the queue on a pod (`examples/runpod_jobs_pod.sh`);
`kink_summary.py --no-tilts`; `candidate_frequency_map.py --screen <the seven screen files, comma list> --kinks
docs/predictions/kink_formation_belt4_polished.json --steps docs/predictions/kink_migration_belt4.json --step-cell 4x2
--step-mover 0 --cheapest --ref vdf9-trfe4-cfe1-14 --out docs/predictions/candidate_frequency_map_belt4_steps_4x2_cheapest.json`;
`pin_design.py`.
NEXT, proposed to the PI: (a) the nitrile with a reversed VDF unit (no chlorine; the nitrile's cut per pin is the best
of the units, the reversed unit places the margin and adds its own cut): screen, polish, kink reference, site scan, all
CPU; (b) the CTFE cell's move from its cheapest site (offset 20) and a path to B through an I-like state; (c) Dipole's
E7 when it lands; (d) the other two helix conformers of the designed rows.

**The leaders' kink sites scanned, and E7 asked of Dipole, 2026-11-14 (PI: "Go on both").** Text:
`docs/CANDIDATE_MATERIALS.md` 5.6 ("The leaders' rates by every kink site") and 8.4 (profiles);
`docs/E1_E2_PREDICTIONS.md` (E7). STATE: DONE, CPU ONLY (spend unchanged: about $11.60 of the $46).
(1) SITE SCANS (`examples/kink_site_scan.py --row ROW [--mirror-kink]`, `kink_site_scan_belt4.json`; `pin_design.py` table
"sites": the rate from the sum of exp(-E/kT) over a lattice's sites and both hands against the 64/29/7 cell's sum):
    64/29/7 (anchor) 364 kHz | 71/29 matrix 33 kHz (ratio 0.092) | AN one per 14: 4.8x, 1.8 MHz, 0.097% at the budget
    CTFE with reversed VDF: 143x, 52 MHz, 0.189% | PVDF with 2 reversed units per 14: 2.9x, 1.1 MHz, 0.081%
The leaders' order stands. AN's cheapest site is 2.15 (mirrored pair, two monomers before the nitrile; the nitrile is a
wall to the other hand, 16.1/8.6). THE CTFE CELL IS CHEAP AT EVERY SITE BUT THE CTFE (0.35 to 3.6): its low kink energy
is the lattice's, and at about kT it says the all-trans state of that cell is barely a lattice: trust it least.
(2) E7 (`sarco/materials/ask-polyfind-dipole-unit-screen-cells-v1.md`, files and README in
`sarco/materials/electrostriction/requests/unit_screen_cells_e7_v1/`; polyfind registration at `7d860b9`): margins of six
cell pairs, single points on the tilted files first, cheapest first:
    a  PVDF, one reversed unit in six `VVVVIV` (72 atoms)      +0.064  (new row `vdf5-rvdf1-6`: 4.1 per reversed unit, as the 14-mer)
    b  VDF/TrFE/AN 70/20/10 (122)  +0.565   b' its 80/20 matrix (120)  +0.229
    c  64/29/7 with AN (170)  +0.197   d  CTFE with reversed VDF (168)  +0.076   e  PVDF, 2 reversed per 14 (168)  +0.162
Prediction file `docs/predictions/e7_polyfind_belt4.json` (`examples/e2_polished_predictions.py --set e7`); read answers
with `examples/crosscheck_e1_e2.py --polyfind "" --dipole "" --e2 docs/predictions/e7_polyfind_belt4.json --e2-dipole
<file>`. It queues behind E2's running relaxation, E5, E5 d and the E2 addendum.
TRAPS. A scan started while another process scans the same row repeats its work (the sites file is locked, the plan is
not): start one process per row and hand. Sites where the kink pair did not stay (`kink_kept` false) are left out of
the sums. The scans are rigid-pendant energies; the one-site tables of 5.6 carry the tilts.
NEXT, proposed to the PI: (a) site scans of the two-unit nitrile cells and of CTFE alone, to see what the reversed unit
adds to CTFE; (b) the single-chain step barriers of the leaders in sixteen chains (a pod-hour), the other rate model;
(c) the other two helix conformers of the designed rows; (d) when E7 a lands, a regio-defect matrix with a loosening
pin in it (PVDF, reversed units, CTFE or AN; no TrFE).

**The unit screen on three axes, 2026-11-14 (PI: "Go on the screen"; scored by the standing objective above).** Text:
`docs/CANDIDATE_MATERIALS.md` 5.6 (and the strain-only notes it put into 3, 5.3, 5.4, 5.5; section 7 says what the 364 kHz
anchor is). STATE: DONE, CPU ONLY, NO POD (spend unchanged: about $11.60 of the $46).
UNITS (`examples/candidate_screen.py`: letters F = TFE, I = reversed VDF, Y = vinyl fluoride, E = ethylene, X =
1,2-difluoroethylene; built from pvdf_belt's CH2 and CF2 and trfe_belt's CHF, no new parameter, junctions unanchored).
One per fourteen in the 71/29 matrix, pin per unit and kink sheet (matrix 6.85):
    reversed VDF -4.8, 5.68 | ethylene -1.7, 6.62 | TFE -1.5, 5.44 | vinyl fluoride +0.4, 6.56 | DFE +2.9, 6.20
    (AN +2.0, 3.06 | CTFE +3.4, 2.04 | CFE +3.7, 3.79 | CDFE +8.8, 2.65)
THREE LOWER THE MARGIN; NONE LOOSENS AS THE BULKY UNITS DO (the cut goes with the all-trans cell's width).
THREE AXES (`examples/pin_design.py`, tables "three_axes" and "combinations" in `pin_design_belt4.json/.md`): relaxation
frequency = 364 kHz x exp(-(E_sheet - 3.79)/kT) (the sheets' population model); strain and tan delta at 34 kHz at the
bias; AND THE AC STRAIN AT THE HEAT BUDGET (`loss.SwitchingChannel.budget_point`: the field solved from 4 MW/m^3, capped at
10 V/um; ceiling 0.19% from the assumed mechanical loss):
    measured film's channel as calibrated   364 kHz   tan d 0.060   0.066% at the budget
    AN one per 14                           1.2 MHz   0.005         0.085%     (one per 28: 52 kHz, 0.43, 0.029%)
    AN with TFE / with ethylene             1.6 / 1.1 MHz  0.004 / 0.019   0.097% / 0.110%   (margins +0.186 / -0.119)
    CTFE with reversed VDF                  52 MHz    0.001         0.189%     (margin +0.076, kink sheet 0.84)
    PVDF, 2 reversed units per 14, no TrFE  0.9 MHz   0.010         0.076%     (margin +0.162, kink 3.24, dipole 0.44)
    the model's CFE cell                    364 kHz   0.000         0.004%
STRAIN ALONE RANKS THEM BACKWARDS: on this rate model the dense weak pin wins and the sparse pin loses. THE TWO RATE
MODELS NOW DISAGREE ABOUT THE DESIGN (dense against sparse); the bench's Vogel-Fulcher fit on a copolymer film decides.
TRAPS. (1) NO BEST DRIVE FREQUENCY: at the heat budget the stroke rate is flat above f_relax, a sweep's maximum is its
ceiling (the exchange optimizer's 277 against 695 kHz; `tests/test_loss.py`). (2) Two units' pins do NOT add (off by
+1.3, -2.7, +1.7, -1.9 kcal/mol per chain): compute the cell. (3) In the two-unit cells the kink pair (site 0) sits ON
the second unit, not like for like with the single-unit cells; no cell's kink sites are scanned, so every kink energy is
an upper bound on the lattice's cheapest. (4) A reversed VDF unit turns its own dipole over: the TrFE-free matrix pays
29% of the dipole. (5) The search's margins of these rows are 0.15 to 0.6 below the polished ones: never read them.
(6) Only two of the four helix conformers ran per row.
REPRODUCE: `examples/candidate_screen.py --sequences ROW --conformations "TTTT,TGTG',TG'TG" --out <screen file>`;
`polish_margins.py --rows ROW`; `kink_migration.py --row ROW --nx 1 --ny 1 --mover all --reference-only`;
`kink_summary.py`; `candidate_loss_map.py --screen <all seven candidate_screen_belt4*.json>`;
`candidate_frequency_map.py --screen <the same> --kinks docs/predictions/kink_formation_belt4_polished.json --ref
vdf9-trfe4-cfe1-14 --out docs/predictions/candidate_frequency_map_belt4_polished.json`; `pin_design.py`.
NEXT, proposed to the PI: (a) kink-site scans of the leading cells (AN; CTFE with reversed VDF; PVDF with reversed
units) so their rates go by the cheapest site; (b) DFT cells from Dipole, cheapest first: PVDF with reversed units
(VDF's groups only), the AN cell, the CTFE cell with a reversed unit; (c) the single-chain step barriers of the leaders
(a pod-hour) for the other rate model; (d) the other two helix conformers of the designed rows.

**The ranking rewritten, the bias on the kink barriers, and what makes the terpolymer work, 2026-11-14 (PI: "Both", then
mid-turn: "I still feel that there are better molecules out there. What are the aspects of pvdf cfe that make it suited?
How do we optimize those characteristics").** Text: `docs/CANDIDATE_MATERIALS.md` sections 3, 4, 5 (new: 5.1 to 5.5), 6
(nitrile paragraphs), 7, 2.1 (nine more rows), 8.2 and 8.5 (more lattices), 8.3 ("What a bias does to these barriers").
STATE: DONE, CPU ONLY, NO POD (spend unchanged: about $11.60 of the $46).
(1) STATICS ON THE POLISHED MARGINS. Every row is polished now but the two FAN cells (`polished_margins_belt4.json`).
`candidate_loss_map.py` takes `--margins` (default) and writes a margin sweep; `examples/pin_design.py` writes
`pin_design_belt4.json/.md` (pin per unit, units' cost against benefit, uniform spacing, statistical against uniform).
THE WINDOW at 55 V/um is 0.12 kcal/mol per monomer wide: half switched at +0.07, a tenth polar at +0.15. ONLY THE PINLESS
MATRIX IS ON IT (71/29 +0.056: 0.61 polar, 0.18% AC strain; 64/36 +0.083); every cell with a pin at 5-10 mol% is off it
(64/29/7 +0.317, 65/30/5 +0.258, CTFE +0.300, AN +0.197). PIN PER UNIT at one per fourteen: AN 2.0, CTFE 3.4, CFE 3.7,
CDFE 8.8; margin = m0 + s/N. THE OFFSET (m0) IS NOT KNOWN BETTER THAN THE WINDOW'S WIDTH: the model ranks by differences.
(2) UNITS, COST AGAINST BENEFIT (5.2). Kink sheet like for like (`kink_migration.py --row ROW --nx 1 --ny 1 --mover all
--reference-only`, then `kink_summary.py`, which adds the tilts and rewrites `kink_formation_belt4_polished.json`):
matrix 6.85, CFE 3.79, AN 3.06, CDFE 2.65, CTFE 2.04; one CFE per twenty 4.59. CUT PER KCAL OF PIN: AN 1.93, CTFE 1.41,
CFE 0.84, CDFE 0.48. The cut falls as 1/N like the pin (CFE: 43 and 45 over N at 14 and 20), so the ratio is the number
to choose a unit by. The search's "CTFE lattice metastable against kinking (-0.6)" is withdrawn (+2.04 like for like).
(3) THE BIAS (8.3; `examples/kink_field.py dipoles | relaxed | table`, `kink_field_belt4.json/.md`). The mixed cell now
carries a field (`packer.set_field`; `MixedModel.dipole`, `mixedcell.cell_dipole`): dE/dE_applied = -(static + induced
dipole). Activation dipoles of every one-chain move in the three sixteen-chain lattices are under 0.82 e.A: AT 55 V/um NO
BARRIER MOVES BY MORE THAN 0.10 kcal/mol, no rate by more than a fifth. Checked by relaxing and settling under the field
in four chains: within 0.0004 at the bias (PVDF and the terpolymer), second order 0.05 at 550 V/um. The field moves the
balance (0.07 per monomer, the window), not the rate.
TRAPS. The all-trans and the kinked reference cells were polished apart and can be turned half round about the chains
(PVDF's are): `kink_field.py` takes the polar axis from the all-trans cell in the lattice's frame and fixes the sign by
the sheet. A converged metric band's top is read instead of its settled saddle only when the search slid (top more than
0.05 above it): `saddle_of`. At 550 V/um the terpolymer's shallow state I is no longer a minimum: not a first-order
check there. A row being polished has no `margin` yet: `pin_design.py` skips it.
THREE READINGS OF WHY THE FILM SWITCHES AT 7 mol% CFE WHEN THE MODEL CRYSTAL DOES NOT (5.3): the chlorine's pin overstated
(Dipole's E2 CFE cell decides); the film switches in its long pin-free runs (statistical 7 mol%: 22% of monomers in runs
of forty or more, 0.038% AC strain under a local margin law against the uniform sequence's 0.001% and the calibrated
film's 0.080%); or at the interfaces (tier 4).
NEXT, proposed to the PI: (a) THE UNIT SCREEN of 5.5: units built from groups the potential already has (TFE, reversed
VDF, vinyl fluoride, ethylene, 1,2-difluoroethylene), two numbers each (pin, cut), about half an hour of CPU per unit,
no pod; (b) the helix-wall cell (a long run with the pins at its ends): tests the local margin law and the blocky
sequence; (c) AN's and CTFE's single-chain barriers in sixteen chains (a pod-hour); (d) ask Dipole for an AN-bearing
cell beside the CFE one.

**What N does next: it falls into L, 2026-11-14 (PI: "Go" on "what N does next").** Text: `docs/CANDIDATE_MATERIALS.md` 8.3
("What N does next"). Files: `kink_states_belt4.json` (states N2, B2, NN4, NN12, L; bands N-N2, N-NN4, N-NN12, A-L,
N-L), each sixteen-chain row's `to_N.then` (copolymers) and `to_L` (all three) in `kink_migration_belt4.json`, tables in
`kink_summary_belt4.md`, job logs in `docs/predictions/gpu_jobs_a100_20261114_n3/` and `_n4/`. Test:
`tests/test_kink_step.py::test_n_is_a_doorway_to_l_in_the_copolymers`. STATE: DONE; BOTH PODS DELETED AND CONFIRMED (43
and 26 minutes, $2.51; about $11.60 of the $46 used in all).
SIXTEEN CHAINS, kcal/mol against the sheet, bands in the per-chain metric:
    from N                              64/29/7                         71/29
    back to the sheet                   2.58                            2.95
    N2 (the move again)                 +14.88, 7.2 above N             +12.66, 5.6 (band top, NOT converged)
    a stack neighbour in N too          pair +19.9/+20.0 (2.3 MORE than two apart), 12.5 / 11.7 to join
                                                                        pair +18.73 (0.8 more), 13.2 / 14.1
    B2 (sheet's hand at N's site)       +12.40                          +12.59
    to L                                1.67 above N                    1.61 above N
L = THE KINK PAIR THREE BONDS ON (A MONOMER AND A HALF), THE SHEET'S HAND (gauche bonds 3, 5): +4.97 (64/29/7), +3.74
(71/29), +5.42 (PVDF): the lowest one-chain state off the sheet in all three (under PVDF's B, +5.91). Straight from
the sheet: 15.29, 15.58, 11.82. FROM N: 10.49 and 10.56 (settled), UNDER N's way back, so N goes on to L 4.6 and 9.4
times per return. THE COPOLYMERS' ONE-CHAIN MOVE IS SHEET -> N -> L: barrier 11.40 / 11.89 (N's, unchanged, so the
cheaper-move frequency map stands), back 6.43 / 8.15. PVDF: N -> L is 5.9 above N (13.66); its moves are the step
(11.64, B) or straight to L (11.82).
L was found by accident (a dip on the matrix's N-N2 band) and built from the sheet for the other two with
`settle.sh ROW 4 2 kinked_crystal.state L torsions 0:180,2:180,3:92,5:-71`. New tool: `settle --copy-mover-to CHAIN`
(job keyword `copyto`): a translate of the moving chain put in its state. A SEARCH FROM A BAND TOP THAT SITS NEXT TO AN
END STATE STILL SLIDES; a converged band's top in the metric is the number to take then (the rows' records say which
was used: `highest`). ODD-BOND KINK SITES HAD NOT BEEN LOOKED AT FOR ONE CHAIN: every state before L had its gauche
bonds on even bonds.
NEXT: (a) does a neighbour follow a chain that is in L (`copyto` + `between`), and L's own next move: this is the
advance question now; (b) L and N in four and thirty-two chains; (c) the other odd-bond sites (1,3; 5,7) for one chain
in the three lattices: L may not be the only one; (d) the matrix's unconverged bands (N-N2, N-NN4).

**The way to N, and bands in the per-chain metric, 2026-11-14 (PI: "Go" on "write the preconditioned band on CPU, then a
short pod").** Text: `docs/CANDIDATE_MATERIALS.md` 8.3 (the last two sub-points of "One chain alone in the copolymer
lattices"), `docs/PERFORMANCE_REVIEW.md` 17.10. Files: `kink_states_belt4.json` (states M, N, P, Q; bands), each
sixteen-chain row's `to_N` in `kink_migration_belt4.json`, `candidate_frequency_map_belt4_steps_4x2_cheapest`, job
logs in `docs/predictions/gpu_jobs_a100_20261114_n/` and `_n2/`. Tests: `tests/test_mixedcell.py::test_neb_in_a_metric_...`,
`tests/test_kink_step.py::test_the_copolymers_cheaper_move_is_to_n_and_pvdfs_is_the_step`. STATE: DONE; BOTH PODS DELETED
AND CONFIRMED (35 and 11 minutes, $1.71; about $9.10 of the $46 used in all).
THE CODE: `polyfind.mixedcell.neb(metric=T)`: the band in `w` with `y = T w` (`chain_metric`'s T), tangents, springs,
constraint projections and FIRE in `w`, the reported force and `max_step` still in the scaled variables; no metric =
the old path bit for bit. `kink_migration.py --band-metric` (needs `--surface`), `band_row.sh ... metric`,
`kink_wide_cell.py between` (metric unless `--plain`). MEASURED: PVDF four-chain bands from nothing on a host, 108 and
116 iterations with tops 10.268 and 12.477 (saddles 10.2684, 12.4644) against 171 iterations to 10.345 and >1000 to
12.49 in per-variable units. A CONVERGED BAND'S TOP IN THE METRIC IS THE SADDLE TO 0.01. It does NOT choose the route:
the terpolymer's direct I-B band converges in 106 iterations on a saddle of 15.293, still above the way through M (14.49).
THE SCIENCE, sixteen chains, kcal/mol against the sheet:
    lattice   the step: B, ahead / back     to N: N, ahead / back      the cheaper move
    PVDF      +5.91, 11.64 / 5.73           +7.72, 12.39 / 4.66        the step, 11.64 (1.0e4 /s)
    64/29/7   +11.79, 14.49 / 2.70          +8.82, 11.40 / 2.58        to N, 11.40 (1.5e4 /s)
    71/29     +10.70, 15.00 / 4.30          +8.94, 11.89 / 2.95        to N, 11.89 (6.5e3 /s)
N = a kink pair of the other hand two monomers on (gauche bonds 4, 6: G- T G+). THE WAY TO IT IS FROM THE SHEET, NOT
THROUGH THE WIDE KINK (an image of the I-N band relaxed to the sheet; the direct band A-N then). The terpolymer's 11.40
is one saddle by two routes (direct band; in pieces through the N-like Q at +10.33: 10.63, 11.396). The matrix's top
settled (11.893, one hyperplane) while its plateau behind it had not converged. PVDF's I-N band runs through B over the
step's own saddle. SO: PVDF moves its kink pair one monomer at a time; BOTH COPOLYMERS MOVE IT TWO MONOMERS ON AND TURN
IT OVER, 3.1 CHEAPER THAN STEPPING ONTO THE TrFE UNIT. On the cheaper move the three lattices are within 0.5 (11.4,
11.6, 11.9), the terpolymer fastest. FREQUENCY MAP ON IT: matrix 158 kHz (0.44 of the terpolymer; 0.40-0.43 on every
earlier protocol), PVDF 243 kHz (was 44 MHz on the steps).
NEXT: (a) what N does next: its own next move (to a third site), whether a neighbour follows it (the `--beside`
question of 8.3 for N), where the pair is after two such moves -- the rate of ADVANCE is not the rate of the first
move; (b) the matrix's plateau behind its top (continue `A-N`); (c) thirty-two chains for the copolymers' N;
(d) `device_min_sites` for mixed cells. Bands from here on: `metric`.

**One copolymer chain alone, 2026-11-14 (PI: "Go" on "the terpolymers' single-chain step in sixteen chains").**
Text: `docs/CANDIDATE_MATERIALS.md` 8.3 ("One chain alone in the copolymer lattices"), `docs/PERFORMANCE_REVIEW.md` 17.10.
Files: rows `<row>@0 2x1` and `4x2 mover 0 ahead host flexible` for `vdf9-trfe4-cfe1-14` (64/29/7) and `vdf10-trfe4-14`
(71/29) in `kink_migration_belt4.json`; `vdf10@0 4x2 ...["again_on_device"]`; the states M and N and the bands through
them in `docs/predictions/kink_states_belt4.json`; job logs in `docs/predictions/gpu_jobs_a100_20261114/`; maps
`candidate_frequency_map_belt4_steps_2x1` and `_steps_4x2`. Tests:
`tests/test_kink_step.py::test_one_copolymer_chain_alone_...`, `tests/test_mixed_batched.py::test_batched_form_on_a_copolymer_cell_...`.
STATE: DONE AS FAR AS THE CARD WAS UP; POD DELETED AND CONFIRMED (one A100 for 98 minutes, about $3.60, and $0.06 for a
first pod deleted 98 s in; with the earlier ones about $7.40 of the $46 the PI named).
NUMBERS, kcal/mol against the sheet, one chain's kink pair one monomer ahead, neighbours flexible, cell fixed:
    64/29/7  2 chains (a sublattice)  I 6.33  S1 6.43  top 11.48  B 8.98   11.5 / 2.5
             4 chains                 I 9.04  S1 9.08  top 13.31  B 13.00  13.3 / 0.3
             16 chains                I 9.05  S1 9.11  top 14.49  B 11.79  14.5 / 2.7   (by way of M; direct band's top 15.41 NOT a settled saddle)
    71/29    2 chains                 I 6.47  S1 7.05  top 12.02  B 11.25  12.0 / 0.8
             4 chains                 I 9.53  S1 9.58  top 13.68  B 12.31  13.7 / 1.4
             16 chains                I 9.36  S1 9.53  top 15.00  B 10.70  15.0 / 4.3   (the same saddle by the direct band and through M)
    PVDF     2 / 4 / 16 chains        tops 8.91 / 12.46 / 11.64; the 16-chain row from nothing on the card is the stored one to 1e-4
THE POINT: from the two-chain protocol to one chain in sixteen every barrier rises by about 3 (3.0, 3.0, 2.7), so THE
DIFFERENCES BETWEEN THE MATERIALS STAY (matrix 0.5 above the terpolymer, PVDF 2.6-2.9 below) AND THE STEP-BASED FREQUENCY
MAP DOES NOT MOVE: 71/29 at 155 kHz (was 147), PVDF 44 MHz (was 27), the anchor is the terpolymer. Absolute rates: 83 /s
and 35 /s per copolymer chain against the film's 2e6 /s: the film's channel is not this step.
TWO THINGS FOUND ON THE WAY. (1) M: a third state of the moving chain between I and B (gauche bond 4 sharp, its
partner spread over bonds 2-3): 11.20 / 13.09 (64/29/7, 4 / 16 chains), 11.68 / 13.43 (71/29). A band from I to B runs
over its shoulder: the four-chain bands found the dip after 900-1100 iterations; the sixteen-chain terpolymer band
CONVERGED (force 0.1) with its top at 15.41 and `saddle_polish` from it slid away twice. A CONVERGED BAND'S TOP IS NOT A
SADDLE ON THE LOWEST WAY; every number here is an upper bound in that sense. (2) N: a kink pair of the other hand two
monomers on (gauche bonds 4 and 6, G- T G+), +8.82 (64/29/7), +8.94 (71/29), +7.72 (PVDF) in sixteen chains: for the
copolymers the lowest state found for one chain off the sheet, below I. The bands from I to N did NOT converge (tops
12.5 after 600 iterations for the terpolymer, 14.8 and 13.3 after 300 for the matrix and PVDF): bounds. FOR THE
TERPOLYMER THE BOUND IS 2 BELOW THE STEP'S 14.5: its cheapest way out of the sheet may be N, not the step. OPEN.
TOOLS (all `examples/`): `mixed_forms_check.py` (the loops against the batched form on a row's own cell);
`gpu_jobs/band_row.sh ROW NX NY [from CELL [stored]] [tiled]` (a row from nothing; bands continued in rounds until
converged); `kink_wide_cell.py table | settle | between` and jobs `settle.sh`, `between.sh` (a held state relaxed and
named; a band between any two held states, its saddle accepted only if converged, above both ends and with negative
curvature); `step_wide.sh` / `saddle_wide.sh` take `row ROW` first (dry run on a host: 64/29/7 in eight chains, 4a x b,
I +9.50, B +12.47, first saddle +9.58; not stored). The one-off curation of the route barriers
(`by_way_of_M`, the terpolymer 4x2 row's `I-B.search_left_the_path`) was a scratch script: the fields are described in
the row's `barrier_polished_note`.
TRAPS: a cell under 256 charged sites is evaluated on the host whatever the pod has (four PVDF chains: 240): do not
queue one. The job runner reads a queue line as `script args`, no environment. The mirror runs inside a Monitor that
expires after 30 minutes: re-arm it or the local copy stops. `letters_of` calls 116 degrees gauche (the matrix's
sixteen-chain M shows as "2, 4").
NEXT, in order: (a) settle the way from I to N for the terpolymer in sixteen chains (a longer band, or better a band in
the per-chain metric, 17.10) -- it decides whether the terpolymer's rate is the step's or N's; (b) N as a start: its
own step, and whether a neighbour joins it; (c) thirty-two chains for the copolymers (`step_wide.sh row ...`), since
their barrier RISES from four to sixteen chains; (d) lower `device_min_sites` for mixed cells.

**The card at work, 2026-11-13 night into 11-14 (PI: "Well, we have $46 in credits left - use them if appropriate").**
Text: `docs/PERFORMANCE_REVIEW.md` 17.9 and `docs/CANDIDATE_MATERIALS.md` 8.3. Files:
`docs/predictions/gpu_jobs_a100_20261113/` (every job's log, the device check's result files and summary). Rows in
`kink_migration_belt4.json`: `8x4`, `8x8`, `16x4` new; `4x4` replaced by the card's (the host's B kept as `on_host`);
`8x2["again_on_device"]`; `8x4["second_route"]`. STATE: DONE, THE POD DELETED AND CONFIRMED BY A READ (fourth pod, one
A100-SXM4-80GB, 63 minutes at $2.19/h, about $2.30; the four pods together about $3.70; the PI's credit stood at $46
before this one). Tests: `tests/test_kink_step.py::test_the_step_in_sixty_four_chains_and_the_thirty_two_chain_row_computed_twice`.
HOW IT WAS RUN: not a benchmark script but A QUEUE OF JOBS ON THE BRANCH. `examples/runpod_jobs_pod.sh` (the pod's start
command, fetched from raw.githubusercontent.com at a commit) installs once (28 s), then follows the branch and runs the
next line of `BENCH_QUEUE` (`examples/gpu_jobs/queue_20261113.txt`), each line a script under `examples/gpu_jobs/`
(`device_check.sh`, `step_wide.sh`, `saddle_wide.sh`); the next job is decided by pushing a line, a `stop` line ends
it, nothing is ever sent to the pod, and results are read from its proxy URL under a serving token (scratch only: every
32-hex string is scrubbed from a log before it is committed). Pod: REST `POST https://rest.runpod.io/v1/pods`, image
`runpod/pytorch:2.2.0-py3.10-cuda12.1.1-devel-ubuntu22.04`, ports `8000/http`, env `BENCH_TOKEN`, `BENCH_QUEUE`,
`BENCH_LIMIT_S`, `BENCH_IDLE_S`, `BENCH_LINGER_S`. A job that tests must unset the pod's `POLYFIND_EWALD_DEVICE` for
pytest (the first device check failed two tests that assert the default module is NumPy).
THE DEVICE (17.9): (1) THE TENSOR'S SIX PRODUCTS AS RANK-K UPDATES (`backend.gram`, `cupy.cublas.syrk`): one evaluation
52 ms, 0.13 s, 0.57 s, 3.67 s at 16, 32, 64, 128 chains (was 49 ms, 0.14 s, 0.76 s, 6.00 s; the prediction was 0.2 s and
2.5 s off, it is 0.19 s and 2.3 s). (2) SEVERAL PROCESSES ON ONE CARD DO NOT PAY: 18.7, 23.0, 28.5, 28.4 evaluations per
second for one to four processes at 16 chains, 7.4, 7.6, 8.6, 8.6 at 32 (no MPS daemon; MPS not measured). The card's
"half idle" at 16 chains is not there for a second process to take. THE FARM ANSWER CHANGES IN ONE PLACE: COUNT ONE
PROCESS PER CARD. (3) IN PRODUCTION HALF OF A RELAXATION IS NOT THE EVALUATION: a step of `relax_on_surface` costs 0.25 s,
1.07 s, 7.5 s at 32, 64, 128 chains where the evaluation is 0.13, 0.57, 3.67 s. The rest is the method's dense linear
algebra on the host: an SVD of the whole closure-normal matrix (4 nc x 44 nc) twice a step (38 ms, 0.15 s, 0.56 s on
four cores) and dense products with the metric (44 nc x 44 nc). The normals are block-diagonal by chain when the cell
is fixed (a hyperplane adds one dense row, a free cell a few dense columns). NEXT SPEED ITEM, NO CARD NEEDED: do
`retract`, `normals` and the metric products of `MixedModel.relax_on_surface` chain by chain.
THE SCIENCE (8.3), all PVDF, one chain's kink pair one monomer ahead, kcal/mol against the sheet:
    cell   chains  I       saddle A-I  saddle I-B  B       ahead / back
    4x2    16      9.963   10.660      11.637      5.911   11.64 / 5.73   (host, before)
    8x2    32      9.974   10.660      11.473      5.692   11.47 / 5.78   (host; the card again: 5e-9, 3e-8, 2e-8, 5e-5)
    4x4    32      9.900   10.635      11.620      5.712   11.62 / 5.91   (card; B was the host's)
    8x4    64      9.888   10.628      11.533      5.495   11.53 / 6.04   (card; second saddle by two routes, 1e-7)
    16x4   128     -       -           -           5.436
    8x8    128     -       -           -           5.430
B: THE TWO DOUBLINGS ADD (8x4 measured 5.495, from the 32-chain cells 5.493) and each series closes by 3.7 (along a) and
3.1 (across) per doubling at the last step; both doublings 5.37, tails 0.02 + 0.03: B IS ABOUT 5.3 IN THE DILUTE LIMIT
(was "about 5.4"; sixteen chains overstate it by 0.6). THE BARRIER AHEAD DOES NOT DEPEND ON THE CELL: the second saddle
is 11.64, 11.47, 11.62, 11.53, no trend, not additive (0.08). THE 11.2-11.3 EXTRAPOLATED ON 2026-11-13 FROM THE 8x2 CELL
("the saddle follows B down three quarters as far") IS WITHDRAWN: that was that cell. Barrier ahead 11.5-11.6 (1e4 /s);
the way back grows with the room, 6.0 in 64 chains and about 6.2 against B's dilute limit (9e7 /s, was 2e8). The wide
kink and its saddle gain nothing along the stack and 0.06-0.09 and 0.03 per doubling across it.
WHEN THE 8x4 SADDLE CAME OUT ABOVE THE 8x2 ONE I checked three things before believing it: the same chains respond by
the same amounts in both (chain by chain), the 4x4 cell (its saddle falls by 0.02 only), and the other route to 8x4
through 4x4 (same energy to 1e-7).
NOT DONE, PROPOSED TO THE PI (each is minutes of card now): the candidate chemistries' single-chain step in sixteen
chains (the frequency map rests on two-chain cells for the terpolymers, where one sublattice moves: PVDF went from 8.9
to 11.6 between those two protocols) -- needs a four-chain band of the terpolymer first, then the let-in route of
`kink_wide_cell.py` (which starts at 2x1; it has no 1x1 -> 2x1 step); the diagonal neighbour's step; prefactors from
Hessians; I and the saddles in 128 chains (about 35 minutes of card each at 7.5 s a step, less after the blockwise
algebra).

**The whole evaluation on a device, 2026-11-13 evening (PI, on the first A100 run: "Are you actually using the GPU? I'm
seeing 0% utilization", then "So port the work - we're using a GPU, let's actually use the GPU" and "I terminated - do
the port properly").** Text: `docs/PERFORMANCE_REVIEW.md` 17.6 (the first run), 17.7 (this port), 17.8 (measured).
STATE: DONE AND MEASURED ON AN A100 (third pod, commit 4349f40, 6 minutes, $0.23; files `pod3_*` in
`docs/predictions/ewald_device_a100_20261113/`; 72 device tests passed on the card; suite here 854 passed, 14 skipped --
NOTE the suite can hang AT EXIT, after its summary, on a forked worker of `tests/test_lattice_table.py`'s parallel build
(Python 3.13, fork + ProcessPoolExecutor shutdown; nothing of this work): if a piped run seems stuck after seven minutes,
`py-spy dump` shows `_python_exit`, and killing the worker lets the summary through). ONE EVALUATION: 49 ms at 16 chains
(was 511 ms with only the sums on the card; 2.88 s on the pod's 20 Xeon cores, 59x), 0.14 s at 32 (was 1.14 s), 0.76 s at
64 (59 s on eight host threads), 6.0 s at 128. THE HOST'S SHARE: 19%, 10%, 4%, 1%. THE CARD'S OWN UTILISATION FIGURE while
evaluations run back to back: 49%, 77%, 91%, 92% (the loop form: 5%, 10%, -, 56%). Agreement: 8.6e-15 against the host
loops at 16 chains; the tiling check 2e-13 to 7e-13 kcal/mol at every size. MEMORY: 1.4, 5.5, 21.1, 82.5 GB -- 128 chains
FILL an 80 GB card. The PI watched the card's meter through the run: CPU-only at first (the installation, my host
baseline at 16 chains, my loop-form comparisons), then "100% GPU and GPU Mem usage, that's more like it" at 64-128.
THE POD SCRIPT NOW RUNS NO HOST BASELINE AND NO LOOP-FORM COMPARISON UNLESS ASKED (`EWALD_BENCH_HOST_MAX`,
`EWALD_BENCH_LOOP_FORM`): a rented card's meter should show the card's work.
NEXT ON THE DEVICE, in order of what the numbers say (17.8): (a) SEVERAL PROCESSES ON ONE CARD at 16-32 chains, where one
process leaves the card half idle (not measured; decides how many cards a band needs -- the farm question); (b) THE SIX
TENSOR PRODUCTS AS RANK-K UPDATES (`cupy.cublas.syrk`): 88 -> 51 ms each at 64 chains, 0.85 -> 0.43 s at 128, i.e. about
0.2 s of 0.76 s and 2.5 s of 6.0 s -- testable only on a card; (c) fewer device calls inside `Ewald.terms` for the small
cells. THE FARM ANSWER given to the PI: two 80 GB cards cover the present plan; four to six only if bands or Hessians at
64+ chains become routine; measure (a) and (b) on a rented card before buying.
THE FIRST RUN, in one line: on an A100 the three ported steps (tensor and field, solve, gradient) ran 46x and 65x faster
than 13 EPYC cores at 16 and 32 chains and agreed to 8e-15, but the evaluation was only 5.3x and 9.5x faster because 70%
of it (more, counting the exclusion loops inside the dipole step) was host loops, the card idle meanwhile. Files:
`docs/predictions/ewald_device_a100_20261113/`.
WHAT MOVED (all in `src/polyfind/`, each held to what it replaces by `tests/test_mixed_batched.py`):
(1) PAIR SUMS on the array module: `mixedcell._pair_sums` -- the chain-pair blocks chosen on the host
(`_chain_pair_blocks`), every block through one kernel (`_pair_kernel`) at the nearest image along the chain, the
same-site columns as a second stack with the bonded scales gathered by image. (2) EXCLUSION CORRECTIONS from flat lists
of the bonded pairs of all chains (`ewald.exclusion_pairs`, `exclusion_correction_pairs`,
`charge_dipole_exclusion_pairs`, `CrystalPacker._exclusion_pairs`, `_polarize(pairs=)`): 86 ms -> 2.3 ms at 16 chains,
on the host (11,520 entries). (3) VALENCE terms of all chains at once (`ChainValence.energy_and_grad_batch`).
(4) THE CHAINS built in one batch and the gradient carried back to torsions and angles as TORQUES, not through the
Jacobian of the coordinates (`chainjac.oligomer`, `oligomer_vjp`, `aligned_block_vjp`, `backbone_transform_jacobian`;
`MixedModel._build_light`): 112 ms -> 4.4 ms at 16 chains, 432 -> 13 ms at 64 (this container). (5) `MixedModel.closure`
no longer builds every chain's whole Jacobian to read two windows of it: 112 ms -> 10 ms per call at 16 chains (432 -> 51 ms
at 64), BIT FOR BIT the same numbers (and 4.6 / 14 ms beside the batched evaluation, to rounding). THIS ONE SPEEDS UP
EVERY CONSTRAINED RELAXATION ON THE HOST TOO, with no change in any number. (6) The dipole solve
is CHOLESKY on a device (`backend.cholesky_solve`; a fifth of the LU solve's time on the A100), with the factorization's
status asked for -- CuPy ignores a failed factorization by default -- and the general solve as the fallback.
WHICH FORM RUNS: `MixedModel(xp=None)` follows the Ewald sums: the batched form on the device when they are there
(`POLYFIND_EWALD_DEVICE=cuda`, at least 256 active sites), otherwise THE LOOPS, unchanged and bit for bit what they
were -- every stored number is still the loops'. `POLYFIND_MIXED_BATCHED=1` or `xp=numpy` asks for the batched form on
the host (same numbers to 1e-15); `xp="loops"` forces the loops. `model.timings = {}`
collects the seconds per part.
TO USE IT in a run: `POLYFIND_DEVICE=cpu POLYFIND_EWALD_DEVICE=cuda python examples/kink_migration.py ...`, serial.
RUNPOD (set up 2026-11-13 at the PI's request, guide https://agents.runpod.io/): the API key is a network secret of the
session's environment for host `rest.runpod.io` (the proxy adds it; it is never visible here) -- saved on the PI's
DEFAULT environment, because the session could not be moved to the `polyfind-sarco` environment created for it; the
Runpod plugin and its MCP server were installed in that session's container only (the MCP sign-in is the PI's step).
A pod is created with `POST https://rest.runpod.io/v1/pods` (image `runpod/pytorch:2.2.0-py3.10-cuda12.1.1-devel-ubuntu22.04`,
ports `["8000/http"]`, start command `bash -c "curl -sL .../examples/runpod_ewald_pod.sh | bash"`, env `BENCH_TOKEN`,
`BENCH_LIMIT_S`, `BENCH_LINGER_S`), its results read at `https://<pod>-8000.proxy.runpod.net/<token>/`, and it is
deleted with `DELETE /v1/pods/<id>` and the deletion confirmed by a read. TRAPS: (a) do NOT send `allowedCudaVersions`:
the API's list stops at 13.0, the A100 hosts report 13.2, and the create then fails as "no instances available" for
every card; (b) a container sees 128-256 CPUs and is allotted 13-32: BLAS at the visible count is slower than one
thread (the benchmark reads the allotment); (c) the image's CUDA 12 libraries do not serve the `cupy-cuda13x` wheel the
driver's version selects: the script falls back to `cupy-cuda13x[ctk]` (works; a few hundred MB more to download);
(d) HOST BASELINES ARE CPU WORK WITH THE CARD IDLE AND THE METER RUNNING: the pod script now times the host at the
first size and one thread count only (`EWALD_BENCH_HOST_MAX`, `EWALD_BENCH_HOST_THREADS` for more); (e) a stopped pod
("Exited") is still on the account: delete it.
NOT MOVED, and why: the backbone recursion (a hundred sequential steps: launches would cost more than the steps), the
bonded lists (too small for a device call to pay), `mixedcell.tilt_relaxation`'s own gradient loop (a Python loop over
the tilt variables around each evaluation: the next thing if tilt relaxations of large cells are run on a device), and
the cube itself (17.5: past about a hundred chains the tensor's reciprocal sum, N^2 K, is the wall, and the remedy is a
particle-mesh solve, not hardware).

**The large cells on a device, 2026-11-13 (PI: "Write it an evaluate - I can get a RunPod instance when you're ready",
and "I have ability to order 4-6 A100s for cheap, so I might build a local GPU farm").** STATE: the port is written and
tested here WITHOUT a GPU. NOTHING HAS BEEN RUN ON A DEVICE: there is no device time anywhere in the repository, and
none may be quoted until the pod run comes back. Text: `docs/PERFORMANCE_REVIEW.md` section 17.
NEXT (waits on the PI): on an A100 pod, `bash examples/runpod_ewald_bench.sh` (about half an hour, most of it the
host's own sixty-four-chain evaluations); it writes `ewald_bench_<time>/` with `log.txt` and the JSON files. When they
come back: (a) add 17.6 "Measured on an A100" from the JSON, nothing from memory; (b) the extras decide whether the
six tensor products become `cupy.cublas.syrk` and the solve Cholesky; (c) answer the farm question with those numbers
(what one card buys one relaxation; which cell sizes fit 40 GB and which need 80; whether cards add up, `--all-devices`);
(d) if the share left on the host is what caps the gain, move the pair sums, the chain building and the exclusion
corrections next.
WHAT IS IN: `Ewald.terms` written once against an array module (`ewald.py`; `Ewald(..., xp=)`, `array_module`,
`device_min_sites` 256, `release`, `device_tensor=`); the tensor kept and the dipoles solved on the device in
`CrystalPacker._polarize`; the switch `POLYFIND_EWALD_DEVICE=cuda` (`backend.ewald_backend`, `erfc_of`, `synchronize`);
`tests/device_double.py` (NumPy underneath, CuPy's refusals on top, counts what crosses), `tests/test_ewald_device.py`
(13 tests, run everywhere), four tests in `tests/test_gpu.py` (real device, skipped here);
`examples/ewald_device_benchmark.py` (`--device none|double|auto|N`, `--chains`, `--host-max-chains`, `--host-threads`,
`--all-devices`; the host is timed at one BLAS thread and at every core and a device is compared with the faster) and
`examples/runpod_ewald_bench.sh` (`EWALD_BENCH_DRY=1` runs the whole script without a GPU on the stand-in; done here
on a copy of the tree, exit 0).
MEASURED HERE (host). One evaluation at ONE BLAS thread, the package's setting: 5.3-6.5 s at sixteen chains, 34.5 s at
thirty-two; at FOUR (`POLYFIND_BLAS_THREADS=4`): 4.4-4.7 s and 24.7 s. A CORRECTION of note (5) of the block below and
of what the PI was told on it ("four BLAS threads change nothing"): that comparison never varied the threads, because
`polyfind/__init__.py` sets the loaded BLAS to one thread at import, and `OPENBLAS_NUM_THREADS=4` does not undo that
in a script that imported NumPy first (every script here; so every large-cell job of this week ran NumPy's BLAS on one
thread whatever its command line said); varied properly, four threads are 1.2-1.4x (the solve halves, the tensor
products gain a quarter to a third, the elementwise part nothing).
A LARGE-CELL JOB THAT HAS THE MACHINE TO ITSELF SHOULD SET `POLYFIND_BLAS_THREADS=4`; jobs side by side should not.
The three ported steps are 86% and 94% of an evaluation (cProfile, one thread: 5.77 of 6.68 s, 31.7 of 33.65 s); what
stays on the host is 0.9 s and 1.9 s (pair sums 0.39, chain building 0.26, exclusion corrections 0.19 at sixteen), so
over this host on four threads a card can give at most 5 times at sixteen chains and 13 at thirty-two until those
move too. The host path is bit for bit the previous commit's: `Ewald.terms` on six
reference cases, and energy and gradient of the two-, four- and sixteen-chain cells. On the stand-in the device path
agrees with the host to 3e-15 of each result's size at sixteen chains; 2.0 MB crosses to the device per evaluation and
69 kB back.
Suite 808 passed, 9 skipped (the GPU ones).
TO USE IT in a run, once measured: `POLYFIND_DEVICE=cpu POLYFIND_EWALD_DEVICE=cuda python examples/kink_migration.py ...`,
serial (no `--workers`).
TRAPS. (1) On a machine with a GPU `POLYFIND_DEVICE` defaults to the device in FLOAT32 for the packer's batched kernels,
and the mixed cell's pair sums borrow those: `MixedModel` now refuses such a packer and names the remedy
(`POLYFIND_DEVICE=cpu`); the benchmark and the pod script set it. (2) A process that has used the device cannot fork
workers that use it: several cards want one fresh process per card (`--all-devices` does that). (3) Device numbers
agree to rounding, not bit for bit, and not from run to run (CuPy's weighted `bincount` is an atomic add). (4) CuPy
14.2 wants NumPy 2.0-2.5 and has wheels for CUDA 12 and 13 only. (5) Read from CuPy's source, each of which would have
cost a pod session to find: its `einsum` would have run the neighbour distances as one tiny matrix product per pair
(the device path multiplies and sums instead); `a @ a.T` is a general product there, not `syrk`; an array has no
`setflags`. (6) Timing repeats of one geometry reuse the kept phases (0.5 s of 6.5): the benchmark drops them before
each evaluation; a timing that does not is that half second too fast. (7) The energy's own difference check of the
gradient has a floor of 2e-5 at sixteen chains (pairs crossing the real-space cutoff at erfc ~ 1e-8); for cells the
host is not run on the benchmark holds the device to the stored two-chain cell's energy per chain instead (1e-13 on
the host). (8) The pod script has not been run on a pod. Its dry mode covers everything but the CuPy steps: the CuPy
install, the check that CuPy reaches the device and its libraries, and the fallback to NVIDIA's runtime wheels
(`cupy-cuda12x[ctk]`) when it does not, are untested. (9) The benchmark's machine record is what caught the thread
error above (`"num_threads": 1` under a label that said four): read it.

**The second chain's step, 2026-11-13 (PI: "Go" on "the second chain's barrier").** Text: `docs/CANDIDATE_MATERIALS.md`
8.3, the sub-bullet "But a stepped chain is recalled, not joined" and its table; table: the last step table of
`docs/predictions/kink_summary_belt4.md`; data: rows `vdf10@0 4x2 mover 4 ahead host flexible beside 0` and `... mover 12
... beside 0` of `kink_migration_belt4.json`; tests: `test_a_stepped_chain_is_recalled_not_joined` and the `--beside` part
of the driver's end-to-end test. RESULT (sixteen chains, every chain flexible; energies against the START, the state
with chain 0 one monomer on, 5.91 above the sheet). The two stack neighbours differ. Chain 12 (it has the stepped chain
at +a): NO wide-kink minimum; one saddle 10.96 at which both gauche bonds are in flight (the trailing one at 137 deg);
both ahead +2.44; barrier 10.96 to join, 8.52 to part again. Chain 4 (stepped chain at -a): wide kink +9.40, saddles
10.43 and 13.12, both ahead +2.44; 13.12 / 10.68. The first chain's own barrier to come back is 5.73: it comes back
6500 times for each time the easier neighbour joins (3e4 /s against 2e8 /s at 300 K). A pair parts at 2e6 /s, a hundred
times slower than one chain returns. From the sheet a pair is 8.35 up behind a top of 16.87: 1.5 /s per chain. So in
zero field the sheet does not advance by chains following one another; an advance needs a drive of about 2 kcal/mol
per stepped chain, and what a drive does to these barriers is NOT computed. Not computed either: the diagonal
neighbours (other sublattice; pair energy 9.05), a third chain.
HOW (driver): `--beside 0` takes the stored B of the row `... mover 0 ...` as the start and measures every energy
against it (`row["beside"]`); `--band-from` a row of the SAME cell whose mover is this row's translate now carries the
guide band's whole displacement field (every chain's deviation from the guide's own straight band, translated) onto
this row's straight band: the first band started at a force of 5.7 with a top of 10.7 (84 and 22.7 with the mover's
coordinates alone) and converged in 66 iterations (the single-chain bands took 160). CHECKS: both-ahead is 2.438 from
either neighbour and equals the stored pair state (8.349 above the sheet); chain 4's wide kink 9.3959 from two starts
(one 129, one 198 above; structures within 0.0004 deg); chain 12: the driver's own relaxation from the wide-kink start
ran to the completed step, a descent in steps of at most 0.05 scaled units from the wide-kink geometry does the same
(280 evaluations, energy falling throughout, ending within 0.001 deg of the stored state; `row["slide"]`), and its
band's fixed end `row["I"]` is a point of that slope marked `not_a_minimum`; its saddle, 0.71 below its band's top, was
tested by descents either way (scratch `connect.py`; stored as `polished["connects"]`): one ends in the start, the other
in the completed step. A CONCERTED PASS FOR THE FIRST CHAIN was looked for and not found: chain 12's saddle geometry carried
onto a chain of the flat sheet (scratch `concerted_single.py`) relaxes onto the sequential path and settles on its first
saddle (10.6599; stored as `polished["found_again"]` of the single-chain row): the 11.64 stands.
TRAPS: (1) the BAND TOPS were far off again: chain 12's band stopped at 11.67 (force 0.148) and its saddle is 10.96;
chain 4's second top 13.35 against 13.12. Never quote a band top. (2) One long band over the whole step (17 images,
scratch `one_band.py`) did not converge (forces stuck near 30, images zig-zagging at the ends): two segments, or a
segment to a point on the slope, do. (3) A relaxation that starts in the wide-kink geometry and ends in the completed
step is not a failed relaxation of I: look at the letters and ask whether I exists (scratch `find_wide.py`).

**Minimising on the closure surface, 2026-11-13 (the follow-up of the sixteen-chain step that settled its first
saddle).** Text: `docs/CANDIDATE_MATERIALS.md` 8.3 (the paragraph of that name and its table); code:
`MixedModel.relax_on_surface`, `mixedcell.chain_metric` / `chain_blocks` / `closure_curvature` /
`closure_multipliers` / `tiled_hessian` / `hessian_columns`, `saddle_polish(on_surface=True)`; driver:
`kink_migration.py --surface` (end states I and B, and `--polish`), `--refine` (relax a stored end state further
as a check; the stored state stays, the check is stored beside it as `relaxed_further`); tests:
`tests/test_mixedcell.py` (three new) and `tests/test_kink_step.py` (the `surface` case and the stored sixteen-chain
checks). RESULT: the sixteen-chain FIRST saddle is 10.660 (band top 10.731; 0.0005 left along the path, 0.00007
across; 87 evaluations, 6 minutes). The SECOND saddle redone the same way is 11.6366, the same to four decimals
(151 evaluations against 818), and the stored I and B relaxed further fall by 0.0000 and 0.0001: the barrier
11.6 / 5.7 and the end states stand. Four chains: both saddles reproduced to four decimals (10.2684, 12.4644) in 78
and 203 evaluations against 391 and 1360. The loose SLSQP result for the first saddle (10.669, 0.016 left across) is
kept in the row as `polished["earlier"]`. WIDER CELLS (rows `8x1`, `8x2`, `4x4`; end states only): B is +6.19 in 8a x b (4a x b: +6.47), +5.69 in 8a x 2b
and +5.71 in 4a x 4b (4a x 2b: +5.91); I is +10.05 in 8a x b and +9.97 in 8a x 2b (flat). So the next doubling is
worth 0.22 along the stack and 0.20 across it; the gain falls by a factor of 5-7 per doubling along a and 3 across.
(SUPERSEDED 2026-11-14 by the block "The card at work" above: 8x4 measured 5.495, 128 chains 5.43-5.44, dilute about 5.3.) If the two add, 8a x 4b is 5.49 and the series' tails are 0.04 and 0.11: B IS ABOUT 5.4 IN THE DILUTE LIMIT (an
extrapolation; the sixteen-chain 5.91 overstates it by about 0.5). The thirty-two-chain states were started from
the sixteen-chain ones with the sheet let in where that cell is farthest from the moving chain (scratch scripts
`wide_seed.py`, `wide_seed_b.py`: `--seed-from` takes only a row of the same cell); 93 and 95 evaluations each.
THE SADDLES IN THIRTY-TWO CHAINS (8a x 2b; no band there: `saddle_polish` on the surface from the 4x2 band's top
images and tangents let into the wider cell, scratch `polish32.py`; stored as `row["saddles"]`): first 10.6601
(sixteen: 10.6600), second 11.4725 (sixteen: 11.6366). So in thirty-two chains I +9.97, B +5.69, barrier 11.47
ahead / 5.78 back: the second saddle follows B down three quarters as far, the way back does not change, the first
saddle and the wide kink do not depend on the cell. (WITHDRAWN 2026-11-14: in 4x4 and 8x4 the second saddle is 11.62 and 11.53; the barrier ahead is 11.5-11.6 in every cell.) Continued to B's dilute limit: about 11.2-11.3 ahead, 5.8-5.9
back. In thirty-two chains every minimisation spent half of its evaluations between an energy settled in the fourth
decimal and its stopping test (gtol 1e-4, or 40 evaluations of stall): use `--surface-gtol 3e-4` or a shorter
`--patience` there. Three thirty-two-chain jobs at once exceed the machine's memory (one was killed); two fit.
WHY SLSQP WAS SLOW AND STOPPED SHORT (measured, four- and sixteen-chain cells): (a) it tries straight steps, which
leave the closure surface at second order, and off the surface a chain meets its own image: 10.331 at a closure
residual of 0.017 deg where the point returned to the surface is 10.279 (from 10.290; a step of Newton's length),
so its steps stay short and most evaluations are rejected trial points (2.2-2.8 per iteration in the stored
relaxations of four to sixteen chains; on the surface 1.0-1.2); (b) conditioning: the Hessian along the surface (energy's plus the closure's curvature times
its multipliers; without that term it has negative values at a minimum, -3.1 against +2.0 measured) runs from 0.031
to 421 in the scaled variables, 0.018 to 2.1 in the per-chain metric (blocks of placement AND internal coordinates;
internal alone does nothing, 0.08 to 415). The metric handed to SLSQP does not help (174 against 185): (a) is the
waste. On the surface without the metric: 110.
TRAPS. (1) `relax`, `relax_on_surface` and `saddle_polish` SORT the free variables (`np.unique`); `free_sets` of
the driver returns them unsorted. A vector over the free variables (a hyperplane's normal, a metric) must be in the
sorted order. Three hours of "the metric makes the first step explode" were a normal passed in the unsorted order.
(2) The blocks come from the four-chain tiling of the same sheet (a minute) and are cached as
`<row>_<kink>.chain_blocks.npz` in `--checkpoints`; they are the SHEET's, used for every state of the row
(recomputing the moving chain's at the start saves 3 evaluations of 29 and costs 88). (3) The metric is bolder than
SLSQP from a strained start: in the small test model, from 13 kcal/mol above the minimum, it ended in another minimum
1 kcal/mol up. The driver runs the damped descent first (`presteps=300`) and the conformation letters of the moving
chain are stored with every end state: look at them. (4) `relax_on_surface` does not enforce the angle and cell
windows (`info["outside_windows"]` counts leavers; none so far) and has no trust box. (5) A thirty-two-chain
evaluation is 28 s alone on one thread (the 38 s first noted was two such jobs side by side), 5.8 times the sixteen-chain
4.9 s. PROFILED 2026-11-13 (scratch `profile16.py`): `Ewald.terms`, the dense pair and reciprocal sums, is 84% of it
(79% at sixteen chains); the dense solve for the induced dipoles is 9% (7%); chain building 2%. Four BLAS threads change
nothing (26.7 s against 28.2): that function is elementwise NumPy and is not threaded. [WRONG, CORRECTED LATER THE SAME
DAY, block at the top: both of those runs were at one thread; four threads are 1.2-1.4x.] An earlier version of this note
blamed the cubic solve; that was not measured and is wrong. `ewald.py` is NumPy only (the CuPy backend of
`polyfind.backend` covers the batched packing kernels, not this path): the large cells are one function away from a
GPU or from threads, and neither port is made. [THE DEVICE PORT WAS WRITTEN LATER THE SAME DAY, block at the top; not
yet run on a device.]
(6) The container restarted once in this window: background jobs died, the working tree and the scratch directory
survived. Commit results as they land.

**The single-chain step in sixteen chains, 2026-11-12 (PI: "Work on the single-chain step in 16 chains").** Text:
`docs/CANDIDATE_MATERIALS.md` 8.3 (table, the paragraph on saddles, the bullet "A sheet is stiff ...") and 8.5; tables:
`docs/predictions/kink_summary_belt4.md`; data: `kink_migration_belt4.json`, row `vdf10@0 4x2 mover 0 ahead host
flexible` (and `4x1`, `2x2`, `4x2 mover 0,4`, `4x2 mover 0,1`); tests: `tests/test_kink_step.py`. RESULT, kcal/mol against the sheet,
one chain's kink pair one monomer ahead of the other fifteen's, every chain flexible, cell fixed: wide kink I +9.96,
B +5.91, second saddle 11.64 (settled), first saddle 10.73 (its band's top; SETTLED ON 2026-11-13 AT 10.66, block above): barrier 11.6 ahead, 5.7 back
(1e4 /s and 2e8 /s at 300 K, 3e12 /s attempt). Four chains, same protocol and now settled too: +9.85, +7.88, saddles
10.27 and 12.46: 12.5 / 4.6. So the cramped cell overstated B by 2.0 and the barrier by 0.8; the 6.4 guessed on
2026-11-10 is 5.91. Checks: the perturbed sixteen-chain sheet relaxes back to the tiled cell (2e-9 per chain); B by
two routes (from the rigid-chain-neighbour state and from the sheet with the mover's torsions reset) 5.911 both;
conformations as intended in every chain. SIZE SERIES for B: 2a x b 7.88, 2a x 2b 7.51, 4a x b 6.47, 4a x 2b 5.91
(room along the stack axis a is what counts; I is 9.85-10.02 in all four): sixteen chains is NOT the dilute limit of
B (8a x b: 6.19, block above). HOW IT GIVES: the moving chain slides 0.32 A back along its axis, its six neighbours shift up to 0.19 A along
theirs and under 0.1 A aside, torsions up to 11 deg (the single-chain slide of 8.1 again). PAIR: two a-neighbours
ahead together +8.35 (4.17 each) against 11.8 apart: the second chain costs 2.4; a diagonal pair (chains 0 and 1)
+9.05 (4.53 each): the second costs 3.1. The barrier of the second chain's step is NOT computed: that is the
next number (it decides whether a stepped chain is joined or recalled, 2e8 /s back). [COMPUTED 2026-11-13, top block: recalled.]
**SADDLES ARE NOW SETTLED, NOT READ OFF BANDS** (`MixedModel.relax(linear=)`, `mixedcell.saddle_polish`,
`kink_migration.py --polish`, stored as `row[seg]["polished"]`, barriers `barrier_forward_polished`): minimum in the
hyperplane normal to the band's tangent at its top, the hyperplane moved by secant steps until the force along the
tangent vanishes; under 0.01 along and across, curvature along negative. WHY: a band's force criterion does not
settle a top that is soft across the path. Four-chain second top: 12.66 (force 0.16, the old stored value), 12.49
after 600 more FIRE iterations, 12.46 settled. Sixteen-chain second top: 11.84 at the 0.1 criterion, 11.64 settled.
Every band of the file was polished: PVDF one sublattice 7.54 / 8.91 (band 7.55 / 8.94), 64/29/7 6.43 / 11.48 (the
"11.4 or 11.7" is 11.48), 71/29 7.05 / 12.02, PVDF in register first saddle 7.38 (band 8.23 at force 1.0; second NOT
settled, its band is too rough). The step-based frequency map was regenerated on them
(`candidate_frequency_map_belt4_steps`; needs `--screen` with all five screen files, `--step-mover 0`, `--ref
vdf9-trfe4-cfe1-14`): 71/29 147 kHz (was 121), tan delta 0.21 at 0.18% (was 0.26), PVDF 27 MHz. OPEN THEN, CLOSED ON 2026-11-13 (10.660, block above): the
sixteen-chain FIRST saddle is its band's top: its polish was stopped by hand after 105 minutes in its first
hyperplane (a second copy with a progress log went 10.731 -> 10.690 in its first 70 evaluations, still descending;
slow, not shown to be wrong); it lies a kcal/mol under the second and does not enter the barrier. The second
saddle's three hyperplanes took 364, 153 and 301 evaluations.
**Tools and traps of the day.** (1) SPEED: one energy+gradient of the 16 x 10 cell 9.9 s -> 4.3 s
(`src/polyfind/ewald.py`: the reciprocal part of the dipole tensor as six symmetric products on half of the k vectors,
off-diagonal components from k_a k_b = [(k_a + k_b)^2 - k_a^2 - k_b^2] / 2; real-space images that cannot hold a pair
inside the cutoff skipped (`Ewald.prune_images`); cos/sin(k.r) kept between the two calls of one evaluation). Every
output bit-for-bit except the tensor (1e-13). NUMPY'S BLAS IS FOUR-THREADED HERE: `OPENBLAS_NUM_THREADS=1` for jobs
side by side. (2) `kink_migration.py`: `--only A|I|B|A-I|I-B` (one piece per process, fields merged under the file
lock, file written whole and moved), `--seed-from ROW`, `--band-from ROW` (mover's internal coordinates from another
row's band), `--checkpoints DIR` (lowest point so far; a band's checkpoint is removed when the band is stored),
`--patience N --patience-tol --end-gradient` (see 3), `--mover a,b` (several chains together), `--fmax` tighter than a
stored band continues it, `--climb-after`, `--band-step`, `--polish`. (3) `relax(patience=, patience_tol=)` stops when
the lowest energy on the closure surface stalls, and every relaxation reports `gradient_max`. SLSQP's own test (energy
change and closure residual under 1e-9 between iterations) was not passed 220 and 80 evaluations after the
sixteen-chain end states' energies had stopped changing in the fourth decimal. IT IS NOT ROUNDING (scatter of that
cell's energy 8e-12; an interim commit, 2e35ba2, said it was and 00d1e09 corrected it); the reason was not looked
into. (4) TRIED AND NOT KEPT: units that even out the curvature per kind of variable (1072 evaluations against 429:
single torsions are the stiffest directions, the soft motions are concerted); a limited-memory quasi-Newton step on
the band force (noisy, images bunch). Both findings are in the docstrings (`relax`, `neb`). (5) A pkill pattern that
matches the row's common prefix kills every job of the row, queued ones included: match on `--only X` or the out file.

**The single-chain barrier with the lattice responding, 2026-11-10 (later; PI: "Work on the relaxed single-chain
barrier").** `examples/chain_slide.py` (tests `tests/test_chain_slide.py`; data `docs/predictions/chain_slide_belt4.json`;
tables in `kink_summary_belt4.md`; text in `docs/CANDIDATE_MATERIALS.md` 8.1, the block after the registry table). One
beta-PVDF chain held at the half-monomer slip in an nx x ny tiling, the rest relaxed. kcal/mol per monomer of the slid
chain: rigid 6.80; neighbours turning 6.02; turning and moving aside 6.07 / 3.70 / 3.36 / 3.26 in 4 / 16 / 36 / 64
chains at fixed cell and 2.86 / 3.09 / 3.12 / 3.13 with the lateral cell free; both series extrapolate (1/chains) to
**3.15**, 46% of rigid. The half-slip is the saddle by the chain's mirror plane (zero axial force on every chain, a
minimum in all else, unmoved when the others' axial shifts are released). Stack neighbours move aside 0.20 A, diagonal
ones 0.17 A and turn 5 deg. Flexing adds nothing, tilts 0.01. The curvature at the minimum is unchanged (14.2); with
the neighbours' axial shifts free it is 5.8 (16 chains), 4.9 (36): anti-plane shear, not a barrier. Short repeats
(one or two monomers per chain) reproduce the ten-monomer crystal exactly, which is what makes 64 chains affordable;
a one-monomer repeat has no internal freedom, `flexible` needs two. Copolymers by the affordable protocol (one chain
in four, lateral cell free; PVDF 2.86): 71/29 2.37 (rigid 5.84), 64/29/7 2.04 (rigid 4.56). A four-chain cell cannot
give way laterally (6.07 against 3.15), so every four-chain number of section 8 with one chain out of step with its
neighbours is an upper bound. **E5 d, DISPATCHED by the PI on 2026-11-11 ("Go on the DFT check")**: the cheapest DFT
test of the 3.15, E5 c's 24-atom c/2 cell (`beta_2a_slide_c2.xyz`) relaxed with the lateral cell free and the repeat
held. Registered for exactly that protocol (`examples/e5d_prediction.py --period-two` ->
`docs/predictions/e5d_polyfind_belt4.json`; `docs/E1_E2_PREDICTIONS.md`, E5 d): +0.708 kcal/mol per monomer of the
cell (4.515 mEh per cell) at Dipole's repeat against 1.668 rigid there, a +5.1%, b +4.2%, relaxed over the DFT's own
rigid part c 0.45; symmetric start = perturbed starts, no chain turns, pendants untilted even from tilted starts.
`examples/crosscheck_e5.py` recognises the record by geometry (or `beta_d_relaxed` in a table). The specification
Dipole reads is the README of `sarco/materials/electrostriction/requests/beta_axial_registry_e5_v1/`.
New flags: `kink_migration.py --ends-only --trust-tiled`. NOTE: a full-suite run under load left pytest sleeping on a
forked child after printing its summary (766 passed, 5 skipped); kill it if it lingers (the next run, 767, exited clean).
**What it changes in section 8, and a trap.** (a) One chain's kink a monomer ahead of its sheet, neighbours as rigid
chains: +11.51 in four chains, +9.32 in sixteen (intermediate 17.81, 17.36): the four-chain cell overstates the
completed single-chain step by about a fifth; the flexible-neighbour band in sixteen chains is NOT run (a day of
compute: 1600 sites, 10 s per gradient). (b) The lone kink sits at the bottom of the registry well, where the lattice's
give does not reach (curvature unchanged), so it barely moves: neighbours fixed 25.61 in sixteen chains (25.25 in four);
neighbours as rigid chains 22.25 at the four-chain minimum's kink shape (23.07 there).
(c) HOW WELL THE LONE KINK IS HELD, MEASURED (`examples/kink_torsion_scan.py`: the two gauche bonds held at +/-tau, the
rest of a stage relaxed from its stored state; data in the rows' `torsion_scan`): neighbours fixed, four chains: a real
minimum (25.25) behind a rim of at least 1.6 kcal/mol (26.90 at 116 deg, then another branch); neighbours as rigid
chains, four chains: a lip of 0.02 (23.07 at 99 deg, 23.09 at 102-103, then down); sixteen chains: no lip (100: 22.17;
112: 20.21; 120: 17.03; 130: 12.98). So the "rigid bodies" column of 8.1's table is a shelf, not a state. HISTORY OF
THIS POINT, so it is not repeated: `MixedModel.relax` returned the sixteen-chain kink to all-trans; a damped descent
(`dilute_kink.py --descend-to 0.05 --descent-only`) crawled 25.61 -> 22.01 in 800 evaluations (2.3 h) without stopping,
and interim commits (792cd66, 482ffed) read that as a minimum and the reversion as an optimiser's jump. Both readings
were wrong: the stretch is flat, the relaxation was right. A relaxation's reversion and a slow descent are both weak
evidence about a rim; the torsion scan is the test. Confirmed by descent alone all the same: PVDF four chains with
flexible neighbours reverts (23.07 -> 15.6 -> 14.1 -> 0.2 -> 0.03, never rising); 71/29 with neighbours fixed reverts.
Tools kept from the detour: `relax(trust=...)` / `relax_in_basin` (boxed quasi-Newton rounds, about ten times the
evaluations), `dilute_kink.py --descend-to/--descent-only/--trust/--start-from`. `--start-from STAGE` continues a stored row. Evaluations in a 16 x 10 cell cost 10 s each: budget
accordingly, and seed a scratch `--out` file with the stored all-trans reference or the driver starts a packing search
on every core.

**Correction and results, 2026-11-10 (read before the 2026-11-08 block below; it corrects it).** The PI asked for
the dilute kink and the migration barrier. Both needed a cell whose chains differ: `polyfind.mixedcell` (energy and
analytic gradients equal to the packer's to rounding; each chain closed exactly by equality constraints, `relax` with
a damped pre-descent, `descend`, `project`; `tile_state`; a climbing-image band on the closure surface,
`neb(constraint=..., climb_force=...)`; 29 tests; suite 755 passed, 5 skipped). Drivers: `examples/dilute_kink.py`,
`kink_migration.py`, `kink_site_scan.py`, `polish_margins.py`, `kink_summary.py` (tables in
`docs/predictions/kink_summary_belt4.md`), `export_mixed_cell.py`. What they found, in order of weight:
(1) **The packing search biases every chiral row.** It refines a recognised state pattern (all-trans, TGTG') under
the line group of the pattern, which a copolymer sequence does not have, and holds a kinked chain periodic only by a
penalty. Polished (every torsion free, closure exact) the all-trans cells fall 0.05-0.25 kcal/mol per monomer and the
TG cells 0.15-0.72, so margins move, all but one toward the helix: 64/29/7 +0.317 (was +0.110), 71/29 +0.056
(-0.051), 64/36 +0.083 (-0.058), 65/30/5 twenty-monomer +0.258 (+0.055), CDFE-14 +0.685 (+0.402), CFE+CDFE-14 +0.807
(+0.641), CTFE-14 +0.300 (+0.135), AN-14 +0.197 (-0.154: its sign changes, AN is a weak helix pin at the film's
dilution and section 6's reading of it is withdrawn), ten-monomer 70/30 +0.213 (+0.145), 80/20 +0.229 (+0.177),
60/40 -0.199 (-0.329), 70/20/10 CFE +0.868 (+0.966), 70/20/10 CDFE +0.392 (+0.146); PVDF unchanged. The matrix is
still degenerate within resolution; the CFE pin is +3.7 kcal/mol per unit at one per fourteen and at one per twenty
(margin +0.06 + 3.7/N) and the film's cell is no longer on the steep part of the curve at the bias, while the
measured film is. The recommendation "one CFE per fourteen to sixteen"
was read off the search's margins and is SUSPENDED (`docs/CANDIDATE_MATERIALS.md` 2.1 and the note at section 5);
sections 3-6 and the old loss and frequency maps still carry the search's numbers. Every long cell is polished
(`polished_margins_belt4.json`); NOT polished: the ten-monomer rows not listed above (`polish_margins.py --rows ...`,
resumable) and the FAN cells (their tilts are not treated). (2) **A lone kink pair does not exist.** It shortens its chain by 1 A and
the axial registry holds the chain (one PVDF chain slid rigidly: 6.8 kcal/mol per monomer; a sublattice half-slip
relaxed with the cell free 0.51 per monomer of the cell against the DFT seat's GFN1 0.56): 25/23 kcal/mol in PVDF
with the neighbours fixed/rigid at ten monomers, 21.7 at infinite repeat, 21/20 in 64/29/7, and it reverts to
all-trans when the neighbours are flexible (71/29: even when fixed). Bias, measured: the mixed cell's chains are
twice too stiff along the axis (674 GPa against DFT's 303); with the angle restraint at 100 (327 GPa,
`--angle-stiffness 100 --angle-window 14`) the ten-monomer numbers are 19.8/18.2, so about 16-17 at infinite repeat.
Kink pairs come as sheets (every chain at the same height). The old section 8's "pinned lattice = upper bound for a
dilute kink" is withdrawn. (3) **The sheet, like for like, with tilts:** PVDF 5.06, 71/29 6.85, 64/29/7 3.79
kcal/mol per kink pair (the search's 4.3, 5.5, 1.8 carried the two biases); the chlorine still loosens the matrix by
3.1. One kink pair per fourteen costs +3.8 while the helix is 4.4 below all-trans per chain: kink pairs attract
along the chain; the low-energy objects are runs of helix and trans with walls of about 2 kcal/mol per chain.
(4) **The step** (one chain's kink pair one monomer ahead, two gauche migrations through the wide kink
G+ T T T G-), kcal/mol per moving chain: PVDF one sublattice 7.55 and 8.94 (saddles, converged, reproduced), ending
1.97 up; one PVDF chain in four 10.35 and 12.66, ending 7.88 up (a sheet is stiff; it advances by stacks along the
short axis); every chain together 8.2 per chain = 16.5 per cell (the concerted route costs more than the sequential
one, and in the copolymer cells its band does not settle: not reported). 64/29/7 and 71/29, one sublattice from a
VDF-VDF pair onto VDF-TrFE: first saddle 6.43 and 7.02, second 11.4 (11.4/11.7 in two runs, not converged) and
12.02 (converged), ending 9.0 and 11.3 up, which the sequence-registered model cell overstates. [THESE ARE BAND TOPS;
the settled saddles of 2026-11-12 are in the block at the top: 7.54 / 8.91, 10.27 / 12.46, 7.38 / not settled,
6.43 / 11.48, 7.05 / 12.02, and one chain in sixteen 10.73 (band) / 11.64 ending 5.91 up.] (5) **The sheet sees
the sequence** (`kink_site_scan.py`; five scans complete, both kink handednesses for 64/29/7 and 71/29, one for
71/29 with the TrFE in pairs `VVVVVRRVVVVVRR`; profiles in `kink_summary_belt4.md`; a site's energy is the site's
own: four 71/29 sites re-relaxed from another site's lattice returned to the same energies and cells). In the matrix
a kink pair costs 4.7-5.4 where a VDF unit follows it (PVDF 5.3) and 7.0-7.8 where a TrFE does, 6.5-7.2 with a TrFE in
it, up to 10.4 at a TrFE pair; the mirrored kink reads the same sequence differently (5.1-5.8 and 8.3-8.4). In
64/29/7 every matrix site is 3.4-5.9 (the chlorine's wider lattice removes the TrFE penalty: 3.1 at the same site,
1.2 between the cheapest sites, 3.42 against 4.66), and with a gauche bond on either side of the CFE the sheet costs
11.4-12.1 (G+ T G-) and 8.1-8.3 (G- T G+): neither kink pair sits at a CFE or passes it. Summed over sites,
exp(-E/kT) per repeat and handedness is 1.0e-3 (71/29), 1.1e-3 (71/29 in pairs), 11e-3 (64/29/7). NOT RUN: the
kink pair on the odd bonds; other sequences. (6) **Frequency
map on two rate models** with polished margins, 64/29/7 anchored at 364 kHz
(`candidate_frequency_map_belt4_polished.md`: rate by the sheet's energy; `_steps.md`: rate by the step barrier,
`--steps ... --step-mover 0`): the 71/29 matrix is at 2 kHz on the first by the one site the sheet was put on and
33 kHz by the sum over its sites (9% of the terpolymer's), and 120 kHz on the second (147 kHz on the settled saddles, 2026-11-12); only the matrix is on the
switching curve (polar fraction 0.61), the 64/29/7 cell is pinned off it (+0.32). Which rate model holds
is the bench's Vogel-Fulcher fit on a copolymer film beside the terpolymer's. (7) Asked of Dipole:
`sarco/materials/ask-polyfind-dipole-registry-kink-v1.md`. **E5 and the E2 addendum DISPATCHED by the PI on
2026-11-10 ("Go on E5 and the E2 addendum"); E6 (the PVDF kink sheet) is NOT dispatched.** E5: seven cells built from
Dipole's accepted beta record atom for atom (`docs/predictions/e5_structures/`, copies in
`sarco/materials/electrostriction/requests/beta_axial_registry_e5_v1/`), polyfind's numbers registered twice
(`docs/predictions/e5_polyfind_belt4.json`, `examples/axial_registry.py`): at its own minimum, per monomer of the cell,
sublattice slid c/8..c/2 +0.198, +0.759, +1.455, +1.788; c/2 relaxed with the slip held +0.513 (b +7.7%; +0.498 from a
symmetry-broken start); one chain of four +0.716, +1.699; and in Dipole's cell (atoms within 0.010 A of the DFT's, the
like-for-like set for the single points) +0.172, +0.656, +1.251, +1.533 and +0.671, +1.585. E2 addendum: the polished
70/30 cells and the 64/29/7 pair in the E2 schema with file hashes (`docs/predictions/e2_polished_polyfind_belt4.json`,
copies of the XYZ files in `sarco/.../requests/polished_e2_cells_v1/`). Ingesters, tested: `examples/crosscheck_e5.py
--dipole <result.json ...>` (recognises records by geometry) and `examples/crosscheck_e1_e2.py --e2
docs/predictions/e2_polished_polyfind_belt4.json --e2-dipole <file>` (add `--e2-rigid-pendants` when Dipole's numbers
are single points on the rigid-pendant files: margins +0.055 and +0.140; the `_tilted` files and relaxations go with
+0.213 and +0.317). All registered in `docs/E1_E2_PREDICTIONS.md`
(its notes file is level with it; re-rendering is safe). Dipole at 15:23 UTC: alpha at BFGS 18 (forces and stress
under threshold, energy gate not met), Gamma's 36 representations and dynmat finished, controller validating, not
delivered: WHEN GAMMA LANDS, ingest it against the registered eps_0 prediction first. Next, in order: rewrite sections 3-6 on the polished margin table (and polish the remaining
ten-monomer rows if they are kept); site scans of further sequences at fixed composition, now that the rule is
known (VDF after the kink 5, TrFE after it 7, a TrFE pair 10, a CFE 8-12); partial sheets (8-16 chains) for
nucleation; the wall of a helix run, which is what moves in the relaxor; a step in a cell whose chains do not carry
the same sequence in register.

**Dipole, 2026-11-08 22:33 UTC.** A power outage on Dipole's host killed both DFT successors (alpha at BFGS 37/120
after 27 h, Gamma at representation 22/36); Dipole marked them lost and started fresh ones under the PI's standing
restart instruction, no checkpoint reuse: `pvdf_alpha_ref_fresh_20261108_v2` (ETA uncertain; the v1 pace says about
three days) and `beta_zero_gamma_fresh_20261108_v2` (40-46 h). Hourly lines on sarco main. Nothing from Dipole lands
before 2026-11-10; the stage-8 items that wait on it (task #16) stay parked.

**Correction, 2026-11-08 (read before the candidate paragraph below).** The candidate screen's first pass ran one
of an isotactic sequence's four TG conformers (TGTG'); for every chlorine and nitrile row the other handedness
(TG'TG) is 0.9-1.5 kcal/mol per monomer lower (a multi-start probe confirms it is the model's effect, not the
search's), so the defect rows' margins are +0.57 to +0.97, not near zero, and the defects read as TG-favouring
pins (CFE +7.9, FAN +7.4, CTFE +6.1, AN +3.9 kcal/mol per unit; CDFE does not pin) rather than margin tuners. The
near-degeneracy comes from the TrFE matrix (29-40 mol%), and a fourteen-monomer cell at the film's composition
(64/29/7, one CFE per fourteen) gives +0.11 where the ten-monomer 70/20/10 cell gives +0.97: the model reproduces
the relaxor at the real composition once the cell holds the matrix between the pins. `docs/CANDIDATE_MATERIALS.md`
sections 2-8 are rewritten on this (the 70/20/10 ranking is withdrawn; the kink-pair bracket of section 8 stands
since it does not depend on the margin). `examples/candidate_screen.py` carries the four conformers
(`TG_VARIANTS`, `conformation()` sized to the repeat, `derive()` takes the lowest TG). The 14-mer dynamics: the
71/29 matrix is degenerate (-0.05) and frozen at the drive (kink pair 5.5, channel 0.8 kHz, saturated at the bias);
one CFE per fourteen holds it off saturation (polar fraction 0.24) and cuts the kink to 1.8 (channel at the drive),
giving tan delta 0.07 for 0.125% AC strain at the operating point (`candidate_frequency_map_belt4_anchor14.md`):
the relaxor mechanism in one cell. The 14-mer matrix is degenerate within resolution from 29 to 36 mol% TrFE
(71/29 -0.05, 64/36 -0.06); CDFE at one per fourteen pins harder than CFE (+0.40 against +0.11; pin strengths are
cell properties: CFE +7.9 per unit in the 20 mol% ten-monomer cell, +2.3 in the 29 mol% fourteen-monomer one, CDFE
-0.3 and +6.3) with the fastest lattice (kink 1.1, channel 1.3 MHz) and no strain at the bias; CFE+CDFE per
fourteen adds the two pins (+0.64), so the pin-and-loosener hypothesis fails as posed. The kink pair's handedness
splits 5-9 kcal/mol in the pinned lattice and 0-0.5 in the relaxed one (section 8). The 14-mer brackets with the
pinned-from-own protocol (`--pinned-from-own`; the pinned-from-all-trans start sticks for long repeats): 64/29/7
1.8-3.5, matrix 5.5-7.6; CDFE-14 kink 1.1, CFE+CDFE-14 0.86 (loosest, hardest pinned), AN-14 2.0 with the cell
TT-ground (-0.15: AN favours the all-trans lattice at the film's dilution, FAN about +0.25 with tilts estimated),
so neither nitrile is the relaxor's pin. Deferred: the kink-position scan (the pinned protocol it needs exaggerates
handedness and sticks; the own basin cannot compare positions), the dilute kink (a mixed-conformation cell), the
migration barrier. CTFE at one per fourteen pins as CFE does (+0.135) with less stroke and an all-trans lattice
metastable against kinking (kink -0.6: dynamics undetermined). The pin density is the dial: margin -0.05 + 2.3/N for
one CFE per N in a 30 mol% matrix (one per twenty +0.055, kink 2.3, channel 157 kHz, 0.154% strain at tan delta
0.20; one per fourteen +0.11, kink 1.8, 364 kHz, 0.125% at 0.07), so under the objective the recommended controlled
sequence is a 30 mol% TrFE matrix with one CFE every fourteen to sixteen monomers, strictly periodic
(`docs/CANDIDATE_MATERIALS.md` section 5). All background computations of 2026-11-08 are complete.

**Where things stand, 2026-11-07 (read this first; the detail is in the dated paragraphs below and in `docs/CRYSTAL_STAGE.md`).**
Candidate pass, 2026-11-07 22:50 UTC (`docs/CANDIDATE_MATERIALS.md`, `examples/candidate_screen.py`,
`examples/candidate_loss_map.py`): ten-monomer sequences packed in TTTT and TGTG' under belt4 with pendant
tilts. TrFE sets the TT-TG margin in steps of about 0.3 kcal/mol per monomer per ten mol% (PVDF +0.74,
80/20 +0.18, 70/30 +0.15, 60/40 TT ground); one CFE per ten beside 20 mol% TrFE brings it to +0.03 (the
degeneracy the relaxor needs), CDFE leaves it, CTFE overshoots to TT ground, CFE without TrFE does
nothing; placement of the CFE moves it by 0.03 per monomer (the local-margin spread a controlled sequence
sets by choice). Loss per strain is not discriminated by this pass (tier 3, the barrier per chemistry,
is the gap). Ranking in the document's section 5; chlorine rows are extrapolated chemistry. Section 6 reads
Dipole's slate of other chemistries (VDCN, CNEPO, FANOMe, SCN, VDC, AN ...) against the same screen and against
Dipole's own registry gate and barrier table: the nitrile homopolymers' all-trans crystals are 6-9 kcal/mol per
monomer above their TGTG' forms (not conformational switches), a nitrile defect beside TrFE drives the lattice to a
saturated TT ferroelectric (the opposite of CFE), and Dipole's barrier table is single-bond torsional barriers (on
its loss term every slate candidate scores 2-2.6 times better than the baseline), which Dipole's own reading says do
not decide the loss: the field moves them 1-4%, switching is collective, and the collective number is section 8's. Nothing on the slate
displaces the VDF/TrFE/CFE family for this objective on current evidence (`examples/slate_screen.py`). Section 8
(`examples/kink_formation_screen.py`, `candidate_frequency_map.py`): the kink-pair formation energy in the packed
cell, the first piece of the dynamics the statics could not give, as a bracket between the relaxed kink crystal
(own packing basin) and the pinned all-trans lattice (`refine_crystal(cell_free=False)`): 4.3-9.3 kcal/mol in PVDF,
4.9-7.6 in the TrFE copolymers, 1.2-3.1 in the CFE terpolymer model, 0.4-0.6 with CDFE, so the copolymer-terpolymer
separation (3.7-4.4, a factor 500-1600 in rate) holds at both ends: the terpolymer's channel follows a 34-66 kHz
drive and the copolymer's is frozen there; loss per strain ranks CDFE then CFE at both ends. CTFE (-0.4 to +20.4)
and the nitrile rows (1.3-11.3, 0.4-12.2) are undetermined: free, their lattices relax to a kinked packing; pinned,
they have no room for a kink. Running: the pinned-then-released protocol (`kink_formation_belt4_pinned_released.json`,
the kink crystal in the basin nearest the all-trans cell). Open: the dilute kink (one kinked chain among all-trans
neighbours, a mixed-conformation cell the packer does not have) and the migration barrier.
The mission's two stages are closed and documented: the large-cell response path is 10x faster
(`docs/PERFORMANCE_REVIEW.md` 16), and the crystal terms of the valence potential are fitted to
Dipole's DFT anchors through seven stages. The current preset is `pvdf-dft-valence-belt3`
(fluorine belt sites, the PBE-D3 chain geometry, a fitted backbone-angle restraint; fitted by
bounded least squares, objective 9.11): against Dipole's records beta's cell is -0.1/-0.4/+0.7%,
the relaxed-ion stiffnesses 1.06/1.01/0.91, the Berry polarization 1.00, the polar Berry slope
0.92, the electronic permittivity 0.93-0.95, the alpha-beta gap and the backbone angles within
1.24 tolerances; its E1/E2 sets are `e1/e2_polyfind_belt3.json` (sections 5.18-5.19). The
pre-registered set on `pvdf-dft-valence` still stands as the primary one. Found 2026-11-07 (5.22):
the refinement's rigid CH2/CF2 groups cannot tilt as DFT's do (Dipole's alpha start geometry has
CF2 rock 2.3, wag 2.9, twist 1.3 deg; beta none by symmetry), and under belt3 alpha releases 0.23
kcal/mol per monomer through rigid group tilts alone, so the gap the fit matched carries that much
pendant strain; `polyfind.pendants.tilt_relaxation` evaluates it, `crystal_stage.py --pendant-tilts`
adds it to the gap row, and the refit from the belt3 point is the current preset
`pvdf-dft-valence-belt4` (objective 11.11; the gap back on the anchor with the relaxation inside
the fit, paid by the fluorine well -12%, the H-F pair +8% and the lateral coupling, C_12 now 0.87 of
DFT's and the largest residual at 1.29 tolerances; everything else within 1.1). The stage-8 shape
term starts from belt4. Prediction on record for Dipole's native alpha positions: CF2 wag 2.5, rock
1.1, twist 0.7 deg. Open on our side: the
lateral coupling (C_12 at 0.91, the lateral axis -0.4%) is held by the other anchors and the prior
together and no radial term adds a direction (5.15, 5.18, 5.20); a shape term (hydrogen belt
sites or a CF2 quadrupole) is the next form change and waits for Dipole's rotated-setting cell
(5.17, requested). Open on Dipole's side: the Gamma DFPT and alpha successors are running after a
WSL shutdown killed the first attempts; polyfind's like-for-like prediction for the DFPT tensor is
eps_0 = 2.23 polar, 4.49 lateral, 2.73 chain (5.16, 5.19), posted on the exchange; belt4 gives
2.23, 4.61, 2.73 (5.22), and it does not move when the valence fit's soft C-C stretch constant is
floored at a physical value; alpha's Gamma-route polar eps_0 does (3.47 to 2.93), so the floored
value is the one to quote for alpha (5.22). Suite: 724 passed, 5 skipped (2026-11-07 18:03 UTC). Fitting from
here: `examples/crystal_stage.py --method trf` (one evaluation per variable per iteration, 5.18).

2026-11-04 research note (cloud session, PI-nominated the same day to drive the electrostriction
research): the literature pass is `sarco/materials/computational/papers/2026-11-04_electrostriction_literature_pass.md`
(the giant strain is a field-induced TG/T3G-to-TT switch shared by crystals, OAF and tie molecules;
the relaxor dispersion is interfacial; loss, frequency and strain are one relaxation, with Ea 0.12 eV
and Tf the discriminating number); the plan is `docs/ELECTROSTRICTION_RESEARCH.md` (tiers: intrinsic
M/Q, switching margin and transition strain, switching barrier and tau(T), then the deformable
connectors), and the request to Dipole is `sarco/materials/ask-polyfind-dipole-electrostriction-verification-v1.md`
(E1 DFPT at strained cells, E2 two-phase PBE-D3 cells per chemistry, E3 one NEB, E4 the bench sweep
in Vogel-Fulcher form). New: `polyfind.electrostriction` (intrinsic electrostriction by the strain
and the polarization routes, agreeing to 5e-5 on force-shifted rigid beta; `TwoState` switch,
`transition_strain`: TGTG' to TTTT is +13.4% axial at the ideal angles), `examples/electrostriction.py`,
`tests/test_electrostriction.py`; `pack.polish`, `mechanics.relax` and `_state_at_strain` take
optional `gtol`/`ftol` (defaults unchanged). The PI's reminder that this seat also supports the
microreactor (fluoro-ribosome) programme is answered in section 10 of the plan: `polyfind.chainend`
gives the torsion-gate field per unit type (VDF end: 1.6 V/nm to write gauche, 0.2 to hold trans)
and beta's registry corrugation (0.48 kcal/mol per monomer, 41 pN per monomer peak) for D21/D24;
`examples/chain_end_gates.py`, `tests/test_chainend.py`. Level 1 of the loss model (PI: "go with level
1") is `polyfind.loss`: a bias-tilted Vogel-Fulcher switching channel giving tan delta, the strain per AC
field and the cycle efficiency at any (f, T, E_bias); on published terpolymer inputs the bias-hold point
at 34-66 kHz carries tan delta 0.06-0.12, not the assumed 0.0125, and crosses the 4 MW/m^3 heat budget
near 25 kHz (plan section 11; the bench sweep in VF form replaces every input). Level 2 ("Go") is
`polyfind.kmc`: the units on a lattice with coupling and random margins, rejection-free kinetic Monte Carlo
under a DC plus AC field, barrier and strength calibrated by measurement, small-signal permittivity at any
bias from the lattice's own fluctuations; first results (plan section 12; 20 x 20 and 32 x 32 agree): the bias axis has a knee near
15 V/um on the first parameters, the lattice's apparent activation (0.43 eV) is a fifth of the way from the bare
barrier to the VF slope. The three-state unit (the polar state's dipole points either way; margin `kT ln 2` for a
zero-field polar fraction of one half; the lattice calibrated on the centroid of its own loss peak) removes the
two-state lattice's over-saturation (strength at 55 V/um 45% of zero bias, measured 50%), puts the zero-bias strain
at 2f, and makes the stroke-at-budget surface rise from zero to a maximum near 20 V/um and decline to the bias-hold point
(by a third on the 32 x 32 lattice, six tenths on the 20 x 20) on the first parameters, the position set by the
polar fraction's rise with DC field, i.e. by the bench's quasi-static S(E_dc) curve (E4 item 7); the bias-hold
point is a question for the calibrated lattice. The coupling x spread grid
settles the temperature slope: 0.36-0.37 eV apparent at every disordered point against VF's 1.0, the lattice ordering
into a ferroelectric past J of about 0.03 eV, so the VF slope is an input (`at_temperature(T, tau_target=...)`
recalibrates the barrier at each T to the bench fit) and the lattice's deliverables are the bias and amplitude
dependence at each temperature; the spread is the knob for the suppression (section 12, last block). The PI runs the E4 bench items personally once an LCR meter is in hand (sarco ask, E4 status line).
Needed from the PI: push `kink-switching` and
`vector-repeat-merge`, which exist only on `E:\source\gh\polyfind`.

2026-11-05 update (the mission: a solver an order of magnitude faster than DFT at DFT-like
accuracy, worked in two stages; all on `claude/polymeric-stable-arrangements-uh06b1`).
**Speed** (`docs/PERFORMANCE_REVIEW.md` section 16): the 120-atom copolymer E1 job, which did
not finish in two hours, runs in 359 s on four cores (1,139 s serial), results identical;
neighbour-list Ewald real space and a BLAS reciprocal tensor (1e-13 agreement), one BLAS
thread per process (OpenBLAS's pool made a 360 x 360 solve 200-500x slower), independent
relaxations in forked workers (`mechanics.parallel_map`, `POLYFIND_WORKERS`), no Coulomb
work in the pair kernel under Ewald; serial profile 314 -> 69 s. Left: the dense
Lennard-Jones kernel, the serial reference build. **Accuracy** (`docs/CRYSTAL_STAGE.md`,
preset `pvdf-dft-valence-crystal`, opt-in): the eight nonbonded terms the conformer trainset
left free were fitted to Dipole's QE beta cell, the antipolar alpha cell, the two margins and
the Berry-phase P. Fixed-charge: P 0.115 -> 0.174 (DFT 0.18), alpha's long axis -6.2% ->
-2.7%, beta's polar axis only -4.5% -> -3.7% with the other axes drifting out by about one
tolerance; the isotropic Lennard-Jones form cannot open the polar H...F contact without
opening the lateral F...F one. With induced dipoles (the E1 setting) the polar axis stays at
-4.8% and P overshoots to 0.214, so the pre-registered set on `pvdf-dft-valence` remains the
better E1 prediction and both sets are kept (`docs/predictions/e1_polyfind_crystal.json`,
`e2_polyfind_crystal.json`). Next, in order: fit in the induced-dipole setting (plumb
`polarizable` through `fitting.predict`, or an antipolar branch on the mechanics path) with a
pair-specific H-F term, adding Dipole's six strained cells as anchors when they land.
**Afternoon, 2026-11-05.** Dipole's six +/-1% axial beta cells landed (exchange main 54897c21):
from their stresses the relaxed-ion stiffness is C_11 22.8, C_22 24.5, C_12 5.9, C_33 302.6 GPa
(record `beta_qe_axial_v1`); the pre-registered prediction has C_11 to 1%, C_22 24% low, C_12 47%
low. Stage 2 of the crystal fit (`docs/CRYSTAL_STAGE.md` section 5; the fit path now takes
`polarizable=`, `lj_cutoff=`, a pair-specific H-F Lennard-Jones term through `FFParameters.lj_pairs`,
and the stiffness as anchors) ran in the setting the predictions are made in: P and the margin
inside tolerance, beta's polar axis -5.0% -> -3.9% and no further, lateral +1.5%, C_12 still 57% of
DFT; preset `pvdf-dft-valence-crystal2`, prediction files `e1/e2_polyfind_crystal2.json`, PE's cell
degrades (not an anchor: a warning). The polar-axis scan says why: at the DFT spacing the
electrostatics pulls the axis in by 1.3-2.7 kcal/mol per A and the fits can only inflate the
repulsion. Charge penetration (`EwaldSpec(charge_damping=b)`, exclusion-consistent, tested) and
the PBE chain geometry (C-F 1.386 against the built 1.350) were each probed: -5.0% -> -4.3% and
-4.4% respectively, neither the lever alone; stage 3 (eleven variables, charge penetration free)
ended at the same floor as stage 2 with the damping off (section 5.6): three forms share the
floor, so the next lever is the fluorine's off-site charge reaching the lattice as a site, with the
PBE geometry and Dipole's DFPT tensors as anchors. **E2 correction**: the E2 script omitted the refined
backbone angles' strain energy, favouring the gauche conformers by about 0.45 kcal/mol per
monomer (PVDF TTTT margin +1.30 -> +0.85, copolymer +1.09 -> +0.46); fixed, files regenerated,
Dipole told.
**2026-11-06.** Stage 4 of the crystal fit broke the floor (`docs/CRYSTAL_STAGE.md` 5.7-5.9):
five probes showed that nothing isotropic moves beta's a:b ratio (0.530 against DFT's 0.558) and
that flattening the fluorine's repulsion in the plane perpendicular to C-F does; implemented as a
rigid pendant fragment with two massless, chargeless, dispersion-free "belt" sites 0.7 A either
side of each fluorine (pseudo-element X, `polymers.belt_fluorine`/`pvdf_belt`/`belt_variant`),
a repulsion-only X-H Lennard-Jones wall whose depth is the twelfth variable, inert within its own
chain so the conformer energetics are untouched. Result: objective 65.1 -> 16.3, every cell axis of
beta and alpha within 1% except the two chain axes (+1.0%, +1.4%: the backbone angle the lattice
settles at), P 0.179 (fitted to 0.18; Dipole's Berry 0.193 is the anchor from stage 5 on), margin
-0.760 (-0.726), C_11 25.1 (22.8), C_22 24.6 (24.5), C_12 4.9 (5.9) GPa. Preset
`pvdf-dft-valence-belt`, meaningful only with `--variant belt`. Dipole's Berry repeat landed
(0.19269 C/m^2 at the reference cell; the pre-registered P0 is -27% against it, stage 2 -2%, stage
4 -7%). Next: the belt preset's E1/E2 sets and their cross-check, stage 5 with the 0.193 anchor
from the stage-4 point, the PBE-D3 chain geometry on the belt chain, DFPT tensors as anchors.
**Later, 2026-11-06.** The belt preset's E1/E2 sets are written (`e1/e2_polyfind_belt.json`,
`docs/CRYSTAL_STAGE.md` 5.10; the copolymer E1 takes 1,356 s with the sites): beta against
Dipole cell -0.6/-0.1/+0.8%, C_11 1.10, C_22 1.00, C_12 0.83, Berry P 0.93, eps_inf 0.93-0.94 of
the DFPT electronic tensor; E2 PVDF TTTT margin +0.758 (DFT +0.726). Beta's lateral relaxed-ion
permittivity keeps rising across the stages (2.83 -> 3.25 -> 3.82; the chains' libration) and is
the set's sharpest number for Dipole's DFPT eps_0 to test. Cautions: the copolymer packing lands
in a different minimum (gamma 88 vs 92 deg, dipole tilt 38 vs 23 deg), so its rows compare
minima; PE's shear stiffness is 5x the pre-registered value (C_66 11.5 GPa), so a PE anchor comes
before the belt terms are trusted on a hydrocarbon. Stage 5 (the Berry anchor 0.1927, from the
stage-4 point; `crystal_stage_v5.json`) is running; `examples/crystal_stage.py --geometry pbe`
builds the belt chain at Dipole's PBE-D3 geometry for the stage after. Dipole's three-point polar
Berry pair landed (exchange main 6574b77d): d|P|/d eta_polar = -0.326 C/m^2, the relaxed-ion e_11
on one branch; polyfind's sets are 0.71 (pre-registered), 0.96 (stage 2), 0.87 (belt) of it and
0.98/0.98/0.93 relative to their own P (record `beta_qe_berry_axial_v1`, cross-check lines added;
candidate anchor for the next stage). Stage 6, that anchor on the belt form (`crystalfit.polar_slope`,
`STAGE6`, `--stage 6`; the slope comes from the stiffness block's strained states for nothing and
matches the E1 path to 1e-4), is implemented and queued to start from the stage-5 result
(`crystal_stage_v6.json`); write-up `docs/CRYSTAL_STAGE.md` 5.11.
**Exchange**: the Science 2022 main text makes FA the dehydrochlorinated CFE pair, i.e. the
same `-CH=CF-` unit as DB (the CF3 reading is withdrawn; E2 is six chemistries; record and
devstream amended, sarco a498326c). Dipole at 05:12 UTC: no accepted strained units yet, one
axial relaxation in SCF, alpha relaxing. The E1 report now carries the copolymer's
electrostriction column (intrinsic lattice Q about +3 m^4/C^2 along the polar axis). The first stage-5 run was cut by the session's two-hour background limit at evaluation 292 of 320 (best objective 12.14 from 18.15), and its result file was only written at the end, so the fitter now checkpoints its best-so-far to the output file on every improvement (`crystalfit.fit(checkpoint=...)`, `"complete": false` until the final write; `--x0` resumes from it). Stage 5 landed at objective 11.64 (from 18.15 at the stage-4 point): |P| 0.188 against the Berry 0.193, the alpha-beta gap on the anchor, every residual inside 1.2 tolerances, the chain axes (+0.7%, +1.0%) the one systematic left (`docs/CRYSTAL_STAGE.md` 5.11). Stage 6 (both Berry anchors) landed at 11.19 from 12.81: |P| 0.192 (z -0.10) and the slope -0.299 (z 0.90), the cell and the gap unchanged, the carbon diameter's excess gone (x1.016); not registered, stage 7 decides. Found and closed on the way: the E1/E2 scripts' refinement ran under `refine_crystal`'s default potential (truncated Coulomb, no induced dipoles), leaving c 0.2% short and the angles 0.3 deg closed against the fit path; they now pass the full potential (`refine_kw`) and reproduce the fit's cell to four decimals, and take `--geometry pbe` for the belt chain (5.10, last paragraph). The belt sets predate the fix and are left as they are. The PBE-D3 chain geometry at the same point makes the chain axes longer, not shorter: the lattice opens the backbone angle on either geometry, against the refinement's invented 105 kcal/mol/rad^2 restraint, and with the angles frozen at DFT's values the alpha-beta gap flips sign: the gap anchor was being matched through the alpha chain opening its angles (`docs/CRYSTAL_STAGE.md` 5.12). DFT's angles (beta 114.63 planar; alpha 118.2 at CH2, 116.9 at CF2 from the PBE-D3 start cell) are now anchors and the restraint's stiffness a variable, on the parameter set (`FFParameters.angle_stiffness`, bound for the refinement by `applied()`): stage 7, `--stage 7 --geometry pbe`, ran from the stage-6 result (5.13). It landed at 11.65 from 25.14 at that point: the restraint fitted at 178 kcal/mol/rad^2, alpha's angles on DFT's (118.05, 116.94), beta's 116.0 against 114.6 (its chain axis, +0.75%, is now that angle: the C-C bond is Dipole's), both Berry anchors and the gap held, every cell axis within 0.8%, C_12 (5.0 against 5.9) and beta's lateral axis (-0.6%) the residuals shared with stage 6 (5.13 has the like-for-like table). Registered as `pvdf-dft-valence-belt2` (`fitting.py`; the thirteenth variable is the restraint, bound for every refinement by `applied()`), which means something only on the belt chain at the PBE-D3 geometry (`--variant belt --geometry pbe`); its E1/E2 sets are `e1/e2_polyfind_belt2.json` (5.14): beta on Dipole's records cell -0.3/-0.6/+0.8%, C_11 1.03, C_22 0.99, C_12 0.85, Berry P 1.00, slope 0.92, eps_inf 0.95/0.93/0.94; E2 PVDF TTTT margin +0.726 on DFT's. Beta's lateral relaxed-ion permittivity is 4.39 (the libration: 2.83 -> 3.25 -> 3.82 -> 4.39 across the sets while P is pinned, i.e. the belt forms imply half the isotropic forms' libration stiffness), the number for Dipole's DFPT eps_0 to test. Cautions: the copolymer packing landed in a third arrangement (gamma 92.1, dipoles along b, 85 deg from a; read its rows in the dipole's frame; E2 puts its three conformations within 0.06 kcal/mol), and PE's cell is 9% smaller than under the belt set with eps_inf 3.04 and C_11, C_22 32, 34 GPa: Dipole's PE cell is the anchor needed before the C/H terms are trusted on a hydrocarbon. Suite after the preset: 710 passed, 5 skipped. **2026-11-07:** Dipole's native alpha relaxation and the beta Gamma DFPT job (SCF and D3 done, PH partial, no dynmat) were killed by a WSL shutdown and are terminal FAILED with zero accepted units; Dipole asked the PI for a recovery route and the PI chose fresh successors, relayed on the exchange devstream with the polyfind seat's priority (Gamma eps_0 lateral, alpha angles, PE cell). Stage 8 waits on those; the FA/DB stereogeometry choice for E2 is still open on Dipole's side. Under the PI's 14-hour standing instruction (05:06-19:06 UTC, hourly check-ins as a scheduled task): (1) the belt sites are out of the Ewald sums and the dipole solve (`CrystalPacker._es_idx`; `VDS98["X"] = 0`; belt copolymer energy call 0.64 -> 0.42 s, predictions moved < 4e-6; `docs/PERFORMANCE_REVIEW.md` 16.6); (2) the sensitivity analysis at the stage-7 point (`examples/crystal_stage_sensitivity.py`, `crystal_stage_v7_sensitivity.json`, 5.15): six of twelve directions constrained, no residual addressable alone, and the lateral coupling (C_12, beta b) reachable only by trading the alpha-beta gap to +2.3 tolerances and beta a to +1.2 (verified by a real evaluation, objective 22.0): the form needs a term that couples the lateral axes without moving the gap; (3) the relaxed-ion dielectric tensor by the Gamma route (`examples/gamma_dielectric.py`, 5.16: the fit-path reference, valence terms pinned at it, the backbone angles on the preset's restraint, belt sites as affine combinations of their CF2 group, all-atom Hessian + Born tensors): beta under belt2 eps_0 = 2.24 polar, 4.57 lateral, 2.74 chain against the rigid-chain route's 2.12, 4.39, 2.45; the libration at 46 cm^-1 carries 2.38 of the lateral 2.48; this is the number for Dipole's DFPT, posted on the exchange as a prediction; (4) basins (`examples/copolymer_basins.py`, 5.17): refining every minimum of the packing screen instead of the first finds, on every fitted potential, a second packing of beta with the chains turned 60 deg (dipoles 119 deg from a, gamma 88.5) degenerate with the DFT setting to 0.1 cal/mol per monomer, and the same pair for the copolymer (dipoles 24 and 85 deg from a, degenerate to 0.1 cal/mol), so the copolymer rows of the sets compared arrangements the screen happened to pick; the chain-orientation energy is unconstrained by every anchor in hand, and the DFT energy of beta in the rotated setting was requested on the exchange (sarco dd4b146c; geometries in `docs/predictions/basins_pvdf-dft-valence-belt2_beta_basin{0,1}.xyz`). The belt2 copolymer E1 cell re-run with the inert-site exclusion: 776 s against 1,356 s, numbers within the relaxation noise. `refined_reference(start=, start_params=)` and `e1_predictions.py --basin/--start-params` refine any screen minimum or a given arrangement; the belt2 copolymer in the pre-registered set's arrangement (`e1_polyfind_belt2_tilt24.json`) lines up with that set row by row (5.17), so the sets' copolymer rows compare once the arrangement is named. The crystal fitter has a bounded least-squares method (`--method trf`, residual vector = anchor z rows + ridge components, forward-difference Jacobian, one evaluation per variable per iteration, parallelisable); the stage-7 run from the stage-6 point reached 9.36 in 235 evaluations (31 min; passes Powell's 11.65 at evaluation 144) and its continuation converged at 9.11 in 53 more (`crystal_stage_v7_trf.json`, `crystal_stage_v7_trf2.json`; 5.18): a better minimum of the same objective, every residual within 1.24 tolerances, C_12 z -0.92 (Powell's point: -1.51), gap -0.21, carbon diameter at the valence value, restraint 228; Powell's stage-7 point was local and 5.15's trade-off holds only there. Registered as `pvdf-dft-valence-belt3`; its E1/E2 sets (both copolymer arrangements), Gamma-route eps_0 (2.23/4.49/2.73, within 2% of the posted belt2 prediction) and Jacobian are written (5.18, 5.19): against Dipole's records cell -0.1/-0.4/+0.7%, C_11 1.06, C_22 1.01, C_12 0.91, P 1.00, slope 0.92; E2 PVDF TTTT +0.736; the copolymer's three conformations within 0.06 kcal/mol on both belt presets. At the belt3 point the Jacobian says closing the lateral coupling buys 1.0 in the anchors' sum of squares and costs 0.8 in the other rows (P, beta a, the alpha angles) plus 3.4 of ridge, and a half-weight ridge probe does not move (7.06 -> 7.00, restraint 228 -> 238 only; `crystal_stage_v7_trf_ridge05.json`): the residual is held by the other anchors and the prior together, so a new term is needed, not a weaker prior (5.18); the requested rotated-setting cell is the datum to build it against. Screened against that Jacobian (`examples/stage8_screen.py`, 5.20): the belt distance and diameter and pair-specific F-F, H-H, C-C, C-F, C-H terms all lie in the span of the twelve variables (outside part 0-3% of their length), so no radial term adds a direction; a new term must have a shape (hydrogen belt sites, a CF2 quadrupole) and should wait for the orientation anchor. `--method trf` is the fitter from here. Suite after belt3: 715 passed, 5 skipped.

2026-11-04 update (cloud session, on `claude/polymeric-stable-arrangements-uh06b1`): the five
CPU performance patches of `docs/PERFORMANCE_REVIEW.md` section 15 are applied (commits
31a2456, bdbcb8f, e364e0b, 1700515; section 15.7 has the measurements on the merged tree):
vectorised FFT screen, single-build enumeration, fast chain geometry, cheaper table build,
closed-form chain Jacobian.  Funnel here 61.9 -> 24.9 s (2.5x with the refinement row under
cProfile, 2.1-2.2x without; section 15.7).  Table values move by < 1e-3 kcal/mol (clear any
on-disk `$POLYFIND_TABLE_CACHE`); line-group end points are unchanged to 2e-9; free-
parametrisation end points can stop elsewhere along their flat valleys (objective within
6e-4 kcal/mol per monomer, lattice energy and bend restraint moving up to 1e-2 each in
opposite directions; `refine_crystal(gradient="fd-geometry")` reproduces the previous
route).  The GPU port (G6) stays a patch.  IP: the solver's methods were assessed for
invention disclosure and found to be known art or implementation; the record is
`sarco/shared/ip/polyfind-methods-prior-art-2026-10-16.md`, and one consolidated defensive
row, row 88 of `sarco/shared/ip/idf-candidates.md`, was logged on the PI's instruction.

2026-10-09 update (on `claude/polymeric-stable-arrangements-uh06b1` unless named):
**the polar/antipolar margin is computed** (`docs/POLAR_MARGIN_DFT.md`): on Dipole's own QE
Hamiltonian (PBE-D3(BJ), SSSP precision, 90/360 Ry), fully relaxed, PVDF's gamma-free
antipolar cell is 1.11 meV (0.026 kcal/mol) per monomer BELOW beta, converged to 0.05 meV;
the sign is a +2.2 meV PBE vs -3.4 meV dispersion balance, so "degenerate" is the robust
reading. Our potential agrees to 0.02 kcal/mol. So for PVDF the static margin does not
separate hard from relaxor; barriers do (kink work below). Dipole's CP2K beta geometry is not
a minimum of their QE Hamiltonian (0.40 eV/A, -13.4 kbar, +8.7 meV/monomer); their Berry/Born
columns are evaluated there. **Far-IR** (`docs/PHONONS.md` last section, `polyfind.infrared`):
the "alpha" phonon row was the delta cell (P2_1cn); see the Phonons section. **Tests**: golden
values to 1e-12 relative, digests via stored arrays; 623 pass / 5 skip. **d figures** in
SCREEN.md/BENCHMARK.md rerun under the proper `e`. **Dipole's `physics-vector-repeat`** is
merged on the side branch `vector-repeat-merge` (8217385) only: the complete placed-cell owner
makes `CrystalPacker.energy` 6-11x and `energy_and_grad` 3-4x slower (energies equal to the
ulp), and Dipole's hold stands; asked in the exchange (consumer note 2026-10-09) whether they
accept a batched path held to the owner by a known-answer test. **Switching energetics**
(branch `kink-switching`, d6ce4f6/264115f off 8217385, `docs/SWITCHING.md` there; needs the
merged torsion): a beta kink is compact (180 deg within 1-2 monomers, lowest on the CF2-CH2
bond); a kink row costs 14.4-15.4 kcal/mol (0.63-0.67 eV) per kink, ~60x the resolution but
moved up to ~4 by the stretch choice; a lone kink pair untwists in a 2x2 host (no converged
lone-kink energy); the row hops one monomer over 8.22 kcal/mol (0.357 eV, first-order saddle,
L = 12, fixed cell only); the neutral 180-deg wall is -0.16 mJ/m^2, unresolvable; harmonic TST
with the model's Vineyard prefactor (4.78 THz) gives ~5 MHz per hop at 300 K, an estimate. On
this potential reversal goes by sheets, not single chains. The pinned host itself drifts (C-F
1.400 vs 1.350 built), a known limit carried by every number. That branch also adds an opt-in
`pair_list` to the complete owner, 7-8x faster at 192 atoms, default bit-for-bit unchanged.

2026-09-29 update (all on `claude/polymeric-stable-arrangements-uh06b1`, pushed): the `e`
definition item is closed -- `mechanics.piezoelectric` returns Vanderbilt's proper tensor and
its converse route holds the nominal field (`docs/ELECTROMECHANICS.md` 5.10; Born-fitted
d33 -9.86, d31 +0.42, not the cause of the shortfall). `phonon.relax_all_atom` owns its
acceptance and the screen phonons record is regenerated under it (tables unchanged).
`pack.polish` wraps `dz` inside its objective. `examples/film_d_budget.py`
(ELECTROMECHANICS 5.11): through this crystal's stiffness, d31 = +20 needs a chain-axis `e` of
6.3 C/m^2 (44x |P|), so d31 is a film-compliance property, not a crystal one; d33 = -32 needs
the polar column at -0.50 to -0.76 C/m^2 (2-3x ours), so the discriminating datum is the
provider's relaxed-ion polar column; the clamped-ion transverse differences move d by <1 pC/N.
New `polyfind.film` (ELECTROMECHANICS 5.12, `examples/film_response.py`): exact laminate of
crystal lamellae + passive amorphous layers along the draw axis. Film d31 +2.2..+6.3 (crystal
+0.42) from amorphous incompressibility; raising only the crystal polar column to -0.77
C/m^2 (3.2x) gives film d33 -32, d31 +12..+17.5, d32 +2.5..+3.1 (measured -32/+20/+1.5):
one crystal deficit explains both. Film inputs swept, never fitted; a measured draw modulus
and crystallinity for the same films would make it a prediction. BUT (5.13): the provider's own data
already estimate their crystal polar column at -0.34 (1.42x ours, via their larger Born charges;
kinematics agree, neither crystal relaxes internally under polar strain), which through the
laminate gives film d33 -14, d31 +6..+8: ~43% / 28-41% of measured, all signs right. The rest
(~2.3x) is outside the crystal (amorphous/interface/charge terms the laminate omits).
The d figures in `docs/SCREEN.md`'s response table and `docs/BENCHMARK.md` were rerun under
the proper `e` on 2026-10-08 (old values kept beside them); SCREEN.md's Step 3 correlations,
Rule 2 and candidates table still use the old `e` and are flagged in place.

2026-09-13 evening provider update: an isolated D: worktree now exists at
`D:\sarco-work\polyfind-vector-repeat`, branch `physics-vector-repeat`,
pushed `cb17be5`. It supplies vector-repeat placed-dipole geometry and a
complete six-column affine model record; see its
`deliverables/clamped_ion_vector_repeat_v2/README.md`. Full relaxed-ion
energy/forces/internal strain/C/S remain open. Do not merge source changes
into this C: checkout while the admitted Sarco native jobs depend on its
unchanged source. GPU is reserved; chemistry is CPU-only. The historical
platform/worktree statements below precede this update.

Read this first after a context clear. It is the map; the documents it points to
are the territory. Everything below is committed and pushed on branch
`claude/polymeric-stable-arrangements-uh06b1`; the working tree is clean, no
agents are running, no worktrees exist.

## What this is

`polyfind` searches for stable chain conformations and crystal packings of
semi-crystalline polymers, PVDF and its candidate relatives, and computes their
field and strain response. The novelty is the search: exact dynamic programming
over discrete torsional states instead of stochastic simulation, screw
decomposition to collapse packing to six variables, a tabulated chain-pair
interaction with FFT lattice sums, analytic gradients throughout. `DESIGN.md` is
the architecture; `README.md` is the entry point; `docs/BENCHMARK.md` is the
measured verdict against the goal.

The goal (`/goal`): a polymer solver that fits the field/strain measurement
needs and is 10-100x faster than existing methods. The `/goal` hook is active
and blocks stopping until met; see "Goal status" below for the honest position.

## The collaboration, and the rule that governs it

**The polyfind repository is the communication channel** with "Dipole", the
DFT/MLIP producer working in the sibling repo `E:\source\gh\sarco`
(`materials/gpu_bundle`). Both sides append dated sections to
`docs/REFERENCE_DATA_REQUEST.md`; that file is a live two-way exchange and its
last ~30 sections are the record of the past four days. Read it from the bottom
up before replying to anything. Never write into `sarco` except as an untracked
file; it has staged work in flight. Rebase-conflicts on that document are
routine (both sides append): resolve by keeping both sides, no markers.

Dipole now also edits Polyfind source directly (2270444 added
`chemical_graph.py` and made `Structure.bonds` a derived view). Review their
code diffs, run the suite, and record the verdict in the exchange.

Dipole works on a different machine with ROCm; this machine has no usable GPU
(`docs/PERFORMANCE_REVIEW.md` section 14). Large runs go to RunPod if ever
needed. The user (Niall) also sometimes drives Dipole's subagents directly.

## The governing requirement (Niall, 2026-09-13)

**The sarcomotor's biggest ask is low tan-delta and a high-frequency
response.** Any strain achieved under those conditions makes a workable
polymer; strain magnitude (d33 = -32) is secondary. Consequences: the
piezoelectric-shortfall chase is validation, not the gate; the
polar-versus-antipolar lattice margin (hard ferroelectric vs relaxor) is the
column that matters and it is the one the model cannot resolve (2026-10-09: for PVDF,
PBE-D3(BJ) puts it at -1.1 meV per monomer, degenerate, agreeing with the model, so the
discriminator moves to switching barriers); crystal
phonons (soft modes, lowest optical frequencies) are the intrinsic proxy the
solver can add. The tool computes no loss and no dynamics today. Recorded
in the exchange as the consumer note of 2026-09-13.

## Goal status, in one paragraph

**Speed: met**, well past target (`docs/BENCHMARK.md`): a full response tensor
in 0.3-7.8 s against the baseline's 87-143 s for a single finite-difference
derivative, and the baseline does no structure search at all. **Static side,
shipped preset `pvdf-dft-valence-flux-born` (rechecked 2026-10-08):** clamped-ion
dielectric within 5% of Dipole's DFPT (2.229/2.146/2.520 vs 2.235/2.253/2.601,
polar/long/chain), transverse Born charges within 0.19 e, C33 +6.4%, cell a/b/c
-3.6/+1.7/-1.2%; but |P| 0.143 is 19-24% BELOW the DFT 0.176-0.188 (the earlier
"within a few percent" was never true of this preset). Chain-axis components
remain open in Dipole's convergence ladder; every DFT reference is
`quantitatively_valid=false`.
**Piezoelectric coefficients: short**, crystal d33 = -9.86 vs measured -32, d31 = +0.42
vs +20 (Vanderbilt-proper since 2026-09-29; the earlier -8.8 / +0.02 were the
improper "film" bookkeeping, `docs/ELECTROMECHANICS.md` 5.10). Through the film laminate
(5.12) d31 becomes +2..+6; with the provider's estimated crystal polar column (5.13) the film
is d33 -14, d31 +6..+8, ~40% of measured with every sign right, and the remaining ~2.3x is
non-crystalline. The earlier account of the reason has been narrowed by elimination to one of two things
neither side has: the electronic clamped-ion term, which no classical charge
model represents, or the 0.4-0.6 C/m^2 Berry-slope target itself, which came
from a sweep that failed its own gate. The narrowed request to Dipole is a
**clamped-ion piezoelectric tensor**, which discriminates. Until it arrives, the
tool's honest scope ends at the static side, and that boundary is sharper than
the goal asked for.

## What was eliminated, in order (do not re-run these)

1. Speed was not the constraint; the algorithm was rebuilt around exact DP and
   the pipeline fell from 289 s to 25 s warm. Nine defects found on the way,
   most by reformulating rather than testing (`docs/PERFORMANCE_REVIEW.md`).
2. The potential was unfitted, not merely cheap. A crystal-data fit failed
   (`DESIGN.md` 5.7); a DFT-data fit generalised (5.9); valence terms made the
   forces usable (5.10). The Berry-slope shortfall was never a fitting problem.
3. **Charges**: a Born-consistent charge-flux fit succeeded (4 params, 10
   observations, rms 0.08 e) and moved d33/d31 *away* from measurement, because
   on rigid pendants only the group Born sum acts under strain and it is 0.08 e.
   `docs/ELECTROMECHANICS.md` 5.9.
4. **Polarizability**: literature values, nothing fitted, reproduce the
   transverse dielectric tensor to 3% and close 23% of the d33 gap; cannot reach
   -32 at any physical value. `docs/ELECTROMECHANICS.md` 5.8.
5. **Internal strain / rigid pendants**: Dipole supplied a geometry-only
   Jacobian; the two models differ exactly where predicted (their C-F bonds
   stretch 0.24 A and turn 31 deg per unit strain, ours rigid to 1e-13) and the
   difference carries no polar dipole through either model's charges: +0.33 vs
   +0.34 C/m^2. `docs/INTERNAL_STRAIN.md`. Phrased as a difference between
   models, not real flexibility, at Dipole's request.
6. Electrostatic scale as a single fix: dead, optima 250x apart
   (`docs/ELECTROMECHANICS.md` 5.7).
7. **Rigid chain rotation as a switching / loss proxy (2026-10-02)**: not usable in this cell.
   Clamped-cell in-phase rotation of beta's two chains: barrier 6.05 kcal/mol/monomer at
   40 deg, and the polar state stays a local minimum past 10 GV/m. Cell-relaxed continuation:
   1.5-2.3 kcal/mol/monomer but bound-limited and branch-tracking (not a true MEP), with
   1D-spinodal "coercive fields" of 14-17 GV/m that are just libration stiffness / dipole
   (the phonon column again). Real switching goes by kinks/domain walls along the chain,
   which a one-monomer, two-chain cell cannot hold. `docs/PHONONS.md`, last section.

## Screening verdicts that stand

`docs/SCREEN.md`, consolidated. The model **cannot decide polar vs antipolar**
for most chemistries: gaps sit inside its own 0.27 kcal/mol per monomer
resolution, and Dipole's three tiers (ours, MACE, vertical DFT) disagree on
the VDCN copolymer's sign, all within that resolution. Accessibility of the
polar phase is the one criterion that resolves; by it TrFE leads, which is
what industry uses (a positive control, not a discovery). Nitrile chemistries
(AN, VDCN, FANOME) cannot reach a polar all-trans phase cheaply and are marked
not rankable on polarity. The "nitrile chains misbehave" narrative was
**withdrawn**: it rested on a frozen-angle fit artefact (fixed:
`fit_ris(angles="relaxed")`, `docs/NITRILE_LANDSCAPE.md`) and on Dipole's
finite-chain gates, whose references were never at minima; on requalified
references the chemistries behave alike and the free-azimuth boundary is what
fails. Ewald is implemented and validated (Madelung to every digit); it changed
no polarity verdict, truncation was never the cause.

## Phonons (new, 2026-09-13)

`src/polyfind/phonon.py` gives Gamma-point lattice dynamics on the packer's own
potential: an all-atom energy assembled term for term (known-answer match to
1e-14 against `packer.energy`), analytic gradient, finite-difference Hessian,
acoustic sum rule to 7e-8 without imposing it. `docs/PHONONS.md` has the
results: at the all-atom minimum with stretch minima pinned to built bond
lengths, beta's softest modes are rigid-chain transverse librations at 34 to
49 cm^-1 (~1 THz), alpha 39, gamma 22. No torsional stiffness enters (the
packer's torsion term is a constant), so twist modes are lower bounds; the
constrained references are not all-atom stationary points (fitted stretch
r0 differ from built bond lengths by up to 0.12 A) and show imaginary
rigid-chain modes there. Reading for the requirement: the harmonic Gamma scale is
~1 THz, and that is all it says; it does not establish the absence of an
intrinsic sub-THz limit (Dipole's correction, accepted: no rates, linewidths,
field coupling or switching dynamics are in a Gamma curvature). Pinned-stretch
results are a changed Hamiltonian, labelled separately. **First comparison with a measured
spectrum (2026-10-08, PHONONS.md last section, `polyfind.infrared`):** the "alpha" row above
is the packer's TGTG' reference, which relaxes to P2_1cn (delta-PVDF's polar arrangement),
not alpha's P2_1/c; the model's P2_1/c cell is 0.030 kcal/mol per monomer above it (inside
the resolution) with b 6.4% short. On the P2_1/c cell, the IR-active antiphase libration
(Rabolt and Johnson's 53 cm^-1 rotatory lattice mode, ~60 at 90 K) is 67.5 (68.1 with
Dipole's torsion), +8 high; nothing IR-active lies within 10 cm^-1 of the ~100/102 line;
the brightest mode is an internal Bu at 72.5 / 88.8 whose position is set by the torsion
term. No beta far-IR data below 100 cm^-1 were found. Still open: Dipole's bulk periodic
curvature. `examples/crystal_phonons.py
--polymorphs beta,alpha,gamma` reproduces (gamma takes ~10 min).

## The lattice-curvature column (new, 2026-09-14)

`examples/screen_phonons.py` and `docs/SCREEN_PHONONS.md`: Gamma-point phonons of every
screened chemistry's all-trans polar cell and exact antipolar cell, at gamma = 90 (the
screen's two cells, reproduced to the quoted digits) and with gamma free (the model's own
minima), each re-polished on the phonon packer and relaxed over every atom at fixed cell.
Record in `deliverables/screen_phonons/`. What it says: (1) the softest transverse
rigid-chain mode of the polar cell is 34.5 cm^-1 for beta-PVDF (construction check against
`docs/PHONONS.md`'s 34.0), 22 for CFE, 16 for TrFE-cand at its gamma-free cell, 46 for CDFE,
53-79 for the bulky-pendant candidates; VDCN's polar cell has 13 and 15 cm^-1 nitrile pendant
librations, the softest in the set and not rigid-chain. (2) The nitrile chemistries' rigid
all-trans repeats are 0.3-0.6 A from the all-atom minimum and their "antipolar" cells become
polar (0.24-0.27 e.A per monomer) once the pendants relax: the screen's rigid-chain polarity
column does not survive for AN, VDCN, VClCN-cand. (3) With gamma free, PVDF's +1.021 and
CFE's +1.501 polarity gaps collapse to -0.000 and +0.071 (DESIGN.md section 6 recorded this
for beta only); of the screen's three resolved-polar verdicts only CDFE's survives; after
all-atom relaxation PVDF, CFE and TrFE are degenerate to 0.02 kcal/mol per monomer. The
screen's polarity column should be re-measured with gamma free on both branches (not done;
`antipolar_cell_exact` screens at 90 only). (4) PVDF's gamma-free antipolar competitor
(a rectangular cell with chain 2 above chain 1 along the long axis) is a true minimum, stiffer
than beta (41/49/61 vs 34/40/49 cm^-1) and degenerate with it: beta is not near a soft-mode
instability on this potential, and polar-vs-antipolar for PVDF is a barrier question this
does not answer. Dipole reviewed the first record and found 12 of 36 relaxed cells a few 1e-6
above the requested force tolerance because `relax_all_atom` discarded termination; fixed at
the owner (`phonon.Relaxation`: unrounded force, termination evidence, Newton polish,
`converged`), consumers require acceptance; record regenerated 2026-09-29 (36/36 accepted, 13 via the Newton polish, every tabulated number unchanged to its printed digits). `pack.polish` now wraps `dz` (and the setting
angles) inside its objective, not only afterwards (2026-09-29): the packer's energy is only
periodic in `dz` for about two repeats outside `[0, c)`. The merge hold below is
about merging *their* branch into ours; Dipole's 2026-09-14 notes say their C: checkout is no
longer touched by incoming code (their relay is on D:), so our own source edits on this branch
are safe to push.

## Dipole's Polyfind extension branch (do not merge yet)

`origin/physics-vector-repeat` (e6150ed, 2026-09-13) is Dipole's isolated
extension of our source: a complete placed-cell energy with real Cartesian
torsion (`CrystalPacker.placed_energy_and_grad(P, H, flip)`), general
triclinic images, vector-repeat valence/flux/exclusion gradients, and a
six-column affine clamped-ion record identical to ours on the four shared
columns. Verified here: 640 pass / 5 skip / 0 fail. **Do not merge it into
the shared branch until Dipole releases the dependency**: their C: checkout
tracks our branch and their running native jobs depend on unchanged source.
Read it in a detached worktree. Torsion changes nothing at Gamma for beta
(primitive two-site repeat cannot twist); it raises alpha's torsional modes
by 4-11 cm^-1. Next lattice-dynamics step: beta with a doubled repeat.

## Dipole's open gates (theirs, not ours to chase)

- Berry polarization: **matrix running, three diagonal columns compared.**
  Their reducer is Vanderbilt-proper (`e = dP/de + d_jk P_i - d_ij P_k`); since
  2026-09-29 so is `mechanics.piezoelectric`'s `e` (the numbers below already
  applied the correction in `clamped_ion_columns.py`). Clamped-ion e_y,jj, their frame, both
  P_y < 0, C/m^2, all `quantitatively_valid=false`: xx theirs -0.1675 / ours
  -0.2960; yy -0.1548 / -0.1559; zz +0.1125 / -0.0108. xy shear (x comp.) theirs -0.2905 / ours -0.3272; yz and xz shears
  are not expressible in our cell. Agreement on the polar column and the
  in-plane shear, a 0.12-0.13 gap on both transverse normal strains. Matrix
  is 25/25 complete on their side, still `quantitatively_valid=false` (xz
  amplitude gate fails on a ~1e-8 vector). Ours regenerated
  2026-09-13 by the scratch pattern around
  `examples/internal_strain_jacobian.py::clamped_ion` (the provider record
  file is gone from sarco, so the full example exits 2; call `clamped_ion`
  directly). Our reference P in their frame: charge-only -0.1187, induced
  -0.0248, total -0.1434. Reproducible record with geometry and source
  hashes: `deliverables/clamped_ion_columns/` from
  `examples/clamped_ion_columns.py` (Dipole asked for it; delivered 2026-09-13). **Do not fit to any Berry number until accepted.**
- **Closed 2026-09-29, ours:** the `e` definition. Derived
  (`docs/ELECTROMECHANICS.md` 5.10): an electrode integrates the nominal
  polarization `F^-1 mu / V0` and a voltage fixes the nominal field `F^T E`,
  so the measured coefficient is Vanderbilt's. `mechanics.piezoelectric` now
  returns it as `e` (old ones kept as `e_reference_volume`, `e_improper`), and
  the converse route holds the nominal field, adding the field-stretch stress
  term both routes had omitted. Born-fitted preset: d33 -8.80 -> -9.86,
  d31 +0.02 -> +0.42, routes agree to 0.02 pC/N on both. Not the shortfall's
  cause (3% and 2% of the gaps); the proper dimensional term gives a
  *positive* d31 (-P S_13), so section 5.3's sign argument is withdrawn.
  Examples now label the Vanderbilt d as "film"; recorded deliverables that
  quote the old film numbers (screen, fit_born_flux runs) are stale until rerun.
- Born/dielectric ladder: 4x8x16 passed its raw sum-rule gate (max 0.0065 e);
  the ladder as declared cannot pass because 4x8x8 failed; a new predeclared
  ladder may follow. Transverse eps (2.253, 2.235) is stable; chain-axis
  converging toward ~2.45.
- VDCN copolymer: polar fixed-cell stage completed under PBE-D3 (force 0.0042
  eV/A); antipolar and full-cell stages pending. Our starts (`deliverables/`)
  are what they are relaxing; the staggered variants there did not help and
  should not be spent compute on.
- Held-azimuth finite-chain ladder: size gates still 0/18 even with azimuth
  held. Finite chains do not converge; the periodic route is the only bulk one.
- Transverse stiffness anomaly (their 112 GPa vs our 17 and published ~20):
  still unisolated, queued behind other work; not on our critical path.

## Discipline that has kept this honest

Twelve-plus retractions across both sides, every one recorded in place rather
than overwritten. Rules that earned their keep: compute the like-for-like
quantity before touching a parameter; check a known-answer case first and
suspect the check before the claim when it fails; keep acceptance quantities
(d33, d31) out of every fit; state parameter count against observation count;
phrase results under the producer's declared bounds; never contract provisional
reference data as if validated. When an agent reports a plausible number
traceable to an invented constant, refuse it.

## Housekeeping

- Tests: `export PYTHONPATH="$PWD/src"; python -m pytest tests -q -p no:cacheprovider`, 724 pass, 5 skip (GPU). ~5-7 min on an idle box. Dipole's runtime reports 3 last-digit
  frozen-literal failures in `test_mechanics.py` (PE/alpha/gamma energies at
  1e-15 relative) that do not reproduce here at any commit; machine noise,
  not a defect. Since 2026-10-08 every recorded golden value is held to 1e-12
  relative and the digest tests compare stored arrays (`tests/data/pendant_builds.npz`);
  `docs/PERFORMANCE_REVIEW.md` 15.6 has the per-test record.
- **Never set `POLYFIND_TABLE_CACHE=1`**: it is read as a directory name and once committed 85 MB of caches. It is gitignored now.
- Agents: launch in worktrees; they often pause on their own background jobs and report "waiting"; that is not stuck. Merge, run the suite, push, then remove the worktree.
- Commit attribution: `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>` and the session URL, per the current system reminder.

## If continuing

On our side: the lattice-curvature column exists now (`docs/SCREEN_PHONONS.md`); its next
steps are the gamma-free re-measurement of the screen's polarity column (the alpha/delta
pair, 0.030 kcal/mol apart, is the same unresolved question), and Dipole's bulk periodic
curvature when it lands; the alpha far-IR comparison is done (above). On the piezoelectric side the xx column says the answer is "both":
a real electronic clamped-ion term of ~0.13 C/m^2 and an overstated target.
When the yy column and the accepted matrix land: write the scope boundary
into `DESIGN.md` and `docs/BENCHMARK.md` with the measured size of the
electronic term, and restate the piezoelectric shortfall against the accepted
Berry-derived total rather than the failed sweep's 0.4-0.6. Either way, reply in `docs/REFERENCE_DATA_REQUEST.md`, then re-run
`examples/electromechanics.py --preset flux` if any parameter changes.
