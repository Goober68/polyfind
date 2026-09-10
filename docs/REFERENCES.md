# Sources for the physical claims in this package

Every load-bearing experimental number quoted in `DESIGN.md`,
`docs/CHEMISTRY_EXTENSION.md`, `src/polyfind/pipeline.py`
(`EXPERIMENTAL_CELLS`), `src/polyfind/enumerate.py` (`KNOWN_CHAINS`) and the
comments in `src/polyfind/polymers.py` has an entry here.  Each entry gives the
claim as the package states it, what the literature actually says, the source,
and a verdict:

* **confirmed** - the package's number matches a primary or reputable source.
* **corrected** - the package's number is close but wrong in detail; the prose
  has been fixed, or the entry says what the right number is.
* **wrong** - the claim does not survive contact with the literature.
* **unverified** - no source found that actually supports it.  This is a real
  verdict, not a placeholder; treat these as unsupported until someone checks.

Where sources disagree the entry gives the range rather than picking a winner.

Numbers that would change computed results used to be collected, unapplied, at
the end of this file.  **That batch has now been applied and every table it
touches regenerated**; what each correction did to agreement with experiment,
including the two places where it made agreement worse, is under
[Changes that were re-measured](#changes-that-were-re-measured).  The pending
list is empty.

## Index

| # | Claim | Verdict |
|---|---|---|
| 1.1 | beta-PVDF cell 8.58 x 4.91 x 2.56, rho 1.97 | confirmed |
| 1.2 | alpha/delta-PVDF cell 4.96 x 9.64 x 4.62, rho 1.92 | confirmed |
| 1.3 | gamma-PVDF cell 4.96 x 9.67 x 9.20, rho 1.94 | **corrected and applied** (rho now 1.93; cell is monoclinic, beta ~ 93 deg, which the packing model cannot express) |
| 1.4 | PE cell 7.42 x 4.95 x 2.55, rho 1.00, Pnam | confirmed (a room-temperature determination) |
| 1.5 | PE herringbone, setting angle +/-48 deg from a, chain 2 offset by c/2 | herringbone and angle confirmed (from-a range 41-48.8 deg); the c/2 offset is a parametrisation artifact, **not** the crystal |
| 2.1 | beta = TTTT, alpha/delta = TGTG', gamma/epsilon = T3GT3G' | confirmed |
| 2.2 | beta has equal backbone angles and exactly trans torsions | confirmed |
| 2.3 | the alternately-deflected beta zigzag is a proposal, not a refinement | confirmed |
| 2.4 | alpha's gauche is near +/-45 deg | **unverified** (45 vs 62 deg) |
| 2.5 | the real gamma chain has trans 171-189 and gauche +/-64 deg | **unverified** |
| 3 | beta-PVDF's polarization is "about 0.13 C/m^2 experimentally" | **wrong** (0.13 is a calculation; measured 0.05-0.10) |
| 3.1 | gamma polar but weaker than beta; alpha near zero | confirmed |
| 4.1 | alpha antipolar, beta polar, gamma polar | confirmed |
| 4.2 | delta and epsilon are the polar/antipolar partners of alpha and gamma | confirmed |
| 5 | alpha and beta nearly degenerate, alpha slightly favoured | confirmed |
| 6.1 | PVDC does not adopt a planar zigzag | confirmed (it is a glide TGTG') |
| 6.2 | PVDC's two backbone angles are both 114 deg | **wrong**; 123 and 114 deg, now applied |
| 7.1 | HH/TT defects, a few percent, known to stabilise beta | abundance confirmed, effect **unverified** |
| 7.2 | UFF Lennard-Jones parameters | confirmed |
| 7.3 | PVDF about half crystalline; melt ~450 K | crystallinity confirmed; 450 K **unverified** |
| 7.4 | two-chain cell with chain 2 at (1/2, 1/2) for PE, alpha, beta, gamma | confirmed; an idealisation for gamma |
| 8.1 | isotactic polypropylene chain repeat 6.50 A | confirmed |

Totals over 22 entries: **14 confirmed, 2 corrected, 2 wrong, 4 unverified**.
Entries 7.1 and 7.3 are split verdicts and are counted under the weaker half.

---

## 1. Unit cells and densities

### 1.1 beta-PVDF (Form I) cell

**Claim** (`pipeline.py`, `EXPERIMENTAL_CELLS["pvdf"]["beta (TTTT)"]`, and the
DESIGN 5.2 table): a = 8.58, b = 4.91, c = 2.56 A, density 1.97 g/cm^3.

**Verified**: orthorhombic, space group Cm2m, Z = 2 monomers, a = 8.58 A,
b = 4.91 A, c = 2.56 A.  Ramer & Stiso: *"the lattice constants [a = 8.58 A,
b = 4.91 A, and c = 2.56 A] and the space group Cm2m ... as determined by
Hasegawa et al. were used"*.  Itoh's Table 4.2 quotes the same three numbers as
the reference experimental values.  The crystallographic density computed from
that cell with Z = 2 CH2CF2 units (M = 64.03 g/mol) is 1.972 g/cm^3.

An earlier and slightly different cell exists: Lando et al. give a = 8.47 A,
b = 4.90 A, c = 2.56 A (Ramer & Stiso, Table 1, footnote b).

**Source**: R. Hasegawa, Y. Takahashi, Y. Chatani, H. Tadokoro, *Crystal
Structures of Three Crystalline Forms of Poly(vinylidene fluoride)*, Polymer
Journal **3**, 600-610 (1972), [doi:10.1295/polymj.3.600](https://doi.org/10.1295/polymj.3.600);
quoted verbatim in N. J. Ramer & K. A. Stiso, [arXiv:cond-mat/0508561](https://arxiv.org/abs/cond-mat/0508561)
and in A. Itoh, PhD thesis, [arXiv:1405.5889](https://arxiv.org/abs/1405.5889), Table 4.2.

**Verdict: confirmed.**

### 1.2 alpha-PVDF (Form II) cell, and delta (Form IV)

**Claim** (`EXPERIMENTAL_CELLS["pvdf"]["alpha/delta (TGTG')"]`): a = 4.96,
b = 9.64, c = 4.62 A, density 1.92 g/cm^3.

**Verified**: a = 4.96 A, b = 9.64 A, c (fibre axis) = 4.62 A, beta = 90 deg,
space group P2_1/c, Z = 4 monomers (2 chains).  Itoh's Table 4.2 lists exactly
these as the reference experimental values for Form II, and the same a, b, c for
Form IV (delta, space group P2_1cn) - so pairing alpha and delta under one cell
in `EXPERIMENTAL_CELLS` is right.  Density computed from that cell with Z = 4 is
1.925 g/cm^3.

Note the crystal system is reported two ways.  Because beta = 90 deg exactly, the
cell is pseudo-orthorhombic and several sources call it orthorhombic (e.g.
*"orthorhombic, a = 0.496 nm, b = 0.964 nm, c = 0.462 nm"*, Orientation of PVDF
alpha and gamma crystals in nanolayered films, Table 2); the space group P2_1/c
is monoclinic.  Nothing in `polyfind` depends on the distinction here, because
this cell's beta really is 90 deg; where a reference cell's beta is *not* 90 deg,
as for gamma, the packing model cannot represent it at all (see 1.3).

**Sources**: Hasegawa et al. 1972 (above); Itoh thesis Table 4.2; Ramer & Stiso
[arXiv:cond-mat/0508561](https://arxiv.org/abs/cond-mat/0508561)
(*"The alpha-phase shows a tg+tg- conformation pattern with four formula units
per unit cell and monoclinic space group P2_1/c"*); low-temperature confirmation
a = 4.96, b = 9.64, c = 4.62, all angles 90 deg in
[arXiv:1604.03919](https://arxiv.org/abs/1604.03919), Table 2;
[PMC4379396](https://pmc.ncbi.nlm.nih.gov/articles/PMC4379396/).

**Verdict: confirmed.**

### 1.3 gamma-PVDF (Form III) cell and density

**Claim as it stood** (`EXPERIMENTAL_CELLS["pvdf"]["gamma/epsilon (T3GT3G')"]`,
and the DESIGN 5.2 table): a = 4.96, b = 9.67, c = 9.20 A, density 1.94 g/cm^3.
The density now reads 1.93; see the verdict.

**Verified**: the three edge lengths are right, but two things are missing or
off.

*The cell is monoclinic.*  Itoh's Table 4.2 lists the reference experimental
values for Form III as a = 4.96, b = 9.67, c = 9.20 A, **beta = 93.0 deg**,
space group Cc.  Lovinger's determination gives a = 4.97, b = 9.66, c = 9.18 A,
beta = 92.9 deg.  An earlier orthorhombic cell (Weinhold, Litt & Lando, space
group C2cm) gives a = 4.96, b = 9.58, c = 9.23 A.  So the spread across sources
on b is 9.58-9.67 A and on c is 9.18-9.23 A, and the monoclinic angle of about
93 deg that two of the three carry was not recorded in `EXPERIMENTAL_CELLS`.
An earlier version of this entry called that a reporting omission "not an error
in the search, which fits the cell angle freely".  **That was wrong.**  The angle
the search fits freely is the crystallographic gamma, between a and b; alpha and
beta are 90 deg by construction, because `pack` puts the chain axis along z and
builds a and b perpendicular to it.  gamma-PVDF's unique angle is beta, between a
and the chain axis, so it is precisely the one that cannot be fitted.  The angles
are now recorded in `pipeline.REFERENCE_CELL_ANGLES` together with the fact that
the package approximates this cell as orthorhombic and what that costs (0.14% in
volume; see the re-measurement section at the end).

*The density does not follow from the cell.*  4.96 x 9.67 x 9.20 A with Z = 8
CH2CF2 units gives **1.928** g/cm^3 (1.930 with beta = 93 deg), not 1.94.  The
quoted 1.94 is what the Weinhold orthorhombic cell (4.96 x 9.58 x 9.23) gives,
1.940 - so the tabulated density and the tabulated cell come from different
determinations.

**Sources**: A. J. Lovinger, *Unit cell of the gamma phase of poly(vinylidene
fluoride)*, Macromolecules **14**, 322-325 (1981),
[doi:10.1021/ma50003a018](https://doi.org/10.1021/ma50003a018); values as
quoted in [PMC4379396](https://pmc.ncbi.nlm.nih.gov/articles/PMC4379396/),
Table 4: *"monoclinic, a = 0.497 nm, b = 0.966 nm, c = 0.918 nm,
beta-angle = 92.9 deg"*.  S. Weinhold, M. H. Litt, J. B. Lando, J. Polym. Sci.
Polym. Lett. Ed. **17**, 585 (1979) and Macromolecules **13**, 1178 (1980)
(cited as ref. 7 by Ramer & Stiso for the C2cm alternative).  Itoh thesis
Table 4.2.

**Verdict: corrected, and applied** - edges confirmed; the density now reads 1.93
for the cell as tabulated; and the cell is monoclinic with beta ~ 93 deg, which
`polyfind.pack` **cannot express** (its alpha and beta cell angles are 90 deg by
construction, only the a-to-b gamma being free), so the package approximates the
gamma cell as orthorhombic and `pipeline.REFERENCE_CELL_ANGLES` records both the
real angles and that limitation.  See
[Changes that were re-measured](#changes-that-were-re-measured).

### 1.4 Polyethylene cell, space group and density

**Claim** (`EXPERIMENTAL_CELLS["pe"]`, and the DESIGN 5.2 table): a = 7.42,
b = 4.95, c = 2.55 A, density 1.00 g/cm^3, orthorhombic.  DESIGN 5.1 separately
quotes the chain repeat as c = 2.55 A.

**Verified.**  The structure is orthorhombic with **two chains per cell, each
contributing two CH2 groups (Z = 4 CH2)**, the chains along c, packed in a
herringbone arrangement:

> *"The unit cell contains two chains, each consisting of 2 CH2 groups, giving a
> total of 12 atoms per unit cell.  The all-trans chains extend along the
> crystallographic c-direction (z-axis) and are packed in a 'herringbone'
> arrangement, characterized by the setting angle psi ... alternating from one
> row of chains to another between the values +/-|psi|."*
> ([arXiv:cond-mat/9612092](https://arxiv.org/abs/cond-mat/9612092))

**Space group: #62 in every case, and the symbol tells you the axis setting.**
The mirror is perpendicular to the chain axis, so the position of `m` in the
symbol says which axis the chain runs along:

* **Pnam** - `m` third, so perpendicular to **c**: chain along c, a ~ 7.4,
  b ~ 4.9, c ~ 2.55.  **This is the convention `polyfind` and the polymer
  crystallography use.**  *"The space group is Pnam, with four CH2 groups in the
  unit cell"* ([arXiv:1604.03919](https://arxiv.org/abs/1604.03919)); the same
  symbol is used by Bunn, Busing, Kawaguchi and the ACS Omega review below.
* **Pnma** - the standard setting, `m` second, perpendicular to **b**: chain
  along b.  A tutorial spells the consequence out: *"set Group to 62, the Pnma
  space group for orthorhombic polyethylene ... lattice parameters a, b, c to
  7.48, 2.55, 4.97 ... With this standard crystallographic orientation, the
  polymer chains will be aligned along direction b"* (gamgi.org).
* **Pbnm** is *not* attested for PE in any source read here.  By the same rule it
  would also put the chain along c but with a and b interchanged (a ~ 4.9,
  b ~ 7.4), so "Pbnm with a ~ 7.4" is internally inconsistent.

One outlier to disregard: the vdW-DF paper cited in 1.5 writes *"base-centered
orthorhombic crystal structure (Pna2_1)"*, a non-centrosymmetric group, and the
cell is primitive, not base-centred; treat that as that paper's slip.

**Cell edges.**  a and b contract markedly on cooling, so values scatter with
temperature far more than with method:

| a (A) | b (A) | c (A) | Method / T | Source |
|---|---|---|---|---|
| 7.40 | 4.93 | 2.534 | X-ray, room temperature | Bunn 1939, quoted in [PMC6641976](https://pmc.ncbi.nlm.nih.gov/articles/PMC6641976/): *"an orthorhombic lattice of space group Pnam with a = 7.40, b = 4.93, and c = 2.534 A"*, and again in [doi:10.1371/journal.pone.0006228](https://doi.org/10.1371/journal.pone.0006228) |
| **7.407** | **4.949** | **2.551** | X-ray fibre, room temperature | Busing 1990, quoted in SIBB **7**(13), [doi:10.5821/sibb.v7i13.1622](https://doi.org/10.5821/sibb.v7i13.1622) |
| **7.417** | **4.945** | **2.550** | neutron, HDPE-d4, 300 K | Tashiro 2D-WAND |
| 7.417 | 4.939 | - | neutron, PE-d4, 300 K | Takahashi 1998 |
| 7.42 | 4.96 | - | X-ray, 77 K | Teare 1959 |
| 7.388 | 4.929 | 2.539 | X-ray, 77 K | Kavesh & Schultz 1970 |
| 7.125 | 4.852 | 2.555 | X-ray, 4.5 K | Kawaguchi, Matsui & Kobayashi 1977, Table I |
| 7.121 | 4.851 | 2.548 | neutron, d-PE, 4 K | Avitabile 1975 |

**The package's (7.42, 4.95, 2.55) is a room-temperature determination, not a
composite.**  It rounds Tashiro's 300 K neutron cell (7.417, 4.945, 2.550) to
three digits, and Busing's room-temperature X-ray fibre cell (7.407, 4.949,
2.551) agrees.  Bunn's 1939 c of 2.534 A is the outlier among room-temperature
values, superseded by the later refinements; **c = 2.55 A is right, and an
earlier draft of this file was wrong to call it the top of the range.**

Two things not to cite: Swan 1962 reports thermal-expansion *coefficients*, not
absolute cell dimensions (the paper is closed and its abstract contains no a, b,
c); and the cell 7.64 / 4.480 / 2.543 that circulates as "Bunn's" is Mueller's
1928 C29H60 paraffin cell, which Bunn used only as an indexing starting point.
Macromolecules **23**, 4608 (1990),
[doi:10.1021/ma00223a018](https://doi.org/10.1021/ma00223a018), is Busing's
X-ray study of gel-spun *hydrogenous* fibre - the deuterated work is Avitabile
1975 and Takahashi 1998.

**Density**: *"significantly less than the crystallographic density of 1.00
g/cm3"* and *"taken as 0.997 g/cm3"* (Martin & Passaglia, J. Res. NBS **70A**(3),
221 (1966)).  Computed from the package's cell with Z = 4 CH2 it is 0.995; from
Bunn's, 1.008.  The quoted 1.00 is right to the precision given.

**Verdict: confirmed** - cell, space group, chain count and density all check
out, and the three edges come from one room-temperature determination rather
than being assembled.

### 1.5 The polyethylene herringbone setting angle, and the c/2 offset

**Claim** (DESIGN 5.2): the herringbone arrangement is found *"(setting angles
+/-48 deg from a, chain 2 offset by c/2)"* as the second minimum.

#### The herringbone: confirmed

See the quote in 1.4.  The two chains are related by the operation
`(1/2 + x, 1/2 - y, z)` of Pnam - an a-glide that reverses the sign of the
setting angle.  That sign reversal *is* the herringbone.

#### The axis convention: the trap in this literature

Both conventions are in use and they are **complementary**, theta_a + theta_b =
90 deg:

* from **a**: *"the setting angle theta, defined as the angle between the a (> b)
  axis and the intra-polymer carbon plane"* (Kleis et al.,
  [arXiv:cond-mat/0611498](https://arxiv.org/abs/cond-mat/0611498)).
* from **b**: *"the setting angle is the angle between the projection of planar
  zigzag chains and the b axis"* (Kawaguchi, Matsui & Kobayashi 1977); *"the
  angle between the molecular plane and the b axis"* (ACS Omega 2018,
  [doi:10.1021/acsomega.8b00506](https://doi.org/10.1021/acsomega.8b00506)).

Because the published values cluster near 45 deg, **the overall range looks
almost the same in either convention while individual attributions flip**.  The
same Avitabile 4 K structure appears as 41 deg in Kleis's from-a table and 49 deg
in ACS Omega's from-b table.  This is the single most dangerous thing in this
corner of the literature.

#### The values

| from b (deg) | from a (deg) | Specimen / method | Source |
|---|---|---|---|
| 41.2 | **48.8** | X-ray, room temperature | Bunn 1939 |
| 42.3 | **47.7** | X-ray, 77 K | Teare 1959 (Kleis Table I) |
| 45 | **45** | X-ray, 77 K | Kavesh & Schultz 1970 |
| 44.4-45.1 | ~45 | neutron, 10-300 K | Takahashi |
| **45.5 +/- 3** | 44.5 | melt-crystallized, X-ray 4.5 K | Kawaguchi et al. 1977, Table I |
| **47-48** | 42-43 | solution-grown single crystal | Kawaguchi et al. 1980 |
| **49** | **41** | neutron, 4 K and 90 K | Avitabile 1975 |

(Bold = the value as printed in a source that names its axis; the other column is
the 90 deg complement.)  A modern review gives the envelope: *"The reported
setting angles for orthorhombic PE are within a narrow range of (45 +/- 4 deg)"*
([PMC9795402](https://pmc.ncbi.nlm.nih.gov/articles/PMC9795402/)).

So **measured from a, the published range is about 41 to 48.8 deg**, and
`polyfind`'s 48 deg is close to Bunn's room-temperature value, which is the
warmest determination in the table.  Different authors attribute the spread to
different causes - Kawaguchi et al. to lattice expansion along a and to
morphology (*"Lamellar crystals have a larger setting angle, when their
crystallite size is smaller and when they are thinner and more distorted"*) -
and the trend is not monotonic once both conventions are unified, so treat the
41-48.8 deg spread as the honest uncertainty rather than fitting a law to it.

Not citable: the "44-48 deg" range widely attributed to D. L. Dorset, *The
setting angle of polyethylene - a critique of powder X-ray determinations*,
Polymer **27**, 1349-1352 (1986),
[doi:10.1016/0032-3861(86)90033-9](https://doi.org/10.1016/0032-3861(86)90033-9),
appears only in search-engine summaries; the paper is captcha-walled and was
never read.  Same for P. J. Phillips & H. T. Tseng, Polymer **26**, 650 (1985),
[doi:10.1016/0032-3861(85)90101-6](https://doi.org/10.1016/0032-3861(85)90101-6).

#### The c/2 offset: **not the crystal structure**

An earlier draft of this file said the c/2 offset between the two chains was
"confirmed as the accepted structure", on the strength of a figure caption
(*"Solid (broken) circles represent CH2-units in the plane (a distance c/2 out of
the plane)"*, Kleis et al.).  **That was a misreading and is withdrawn**: the
caption distinguishes the two CH2 groups *within* each chain, which do sit half a
repeat apart, not the two chains from each other.

Worked from the structure instead: Kawaguchi et al. give the carbon site as
`x/a, y/b, z/c = 0.046, 0.065, 0.25` - z fixed at 1/4 by the mirror.  Applying
the Pnam operators to it puts the corner chain's carbons at z = 1/4 and 3/4
**and the centre chain's carbons at z = 1/4 and 3/4 as well**.  The two chains
are at the same two levels; there is **no c/2 stagger**.  The centre chain comes
from `(1/2 + x, 1/2 - y, z)`, a glide with zero c-translation.

This is not a disagreement with what `polyfind` reports, because the two
descriptions are the same structure.  An all-trans chain has its own 2_1 screw:
translating by c/2 and rotating 180 deg about the chain axis leaves it invariant.
So in `polyfind`'s `(setting angle, dz)` coordinates, "chain 2 offset by c/2 at
setting angle phi" and "chain 2 unoffset at setting angle phi + 180 deg" are the
same atoms - and for a planar zigzag those two setting angles describe the same
plane.  The reported c/2 is a coordinate the parametrisation happened to land in,
not a physical stagger.

**Verdict: herringbone and from-a convention confirmed; the 48 deg value is
inside the published from-a range of 41-48.8 deg and close to its
room-temperature end; the "chain 2 offset by c/2" is a parametrisation artifact
rather than a feature of the crystal, and reads as a prediction when it is not.**

**Sources for 1.4 and 1.5**: J. Kleis, B. I. Lundqvist, D. C. Langreth,
E. Schroder, *Towards a working density-functional theory for polymers:
First-principles determination of the polyethylene crystal structure*,
[arXiv:cond-mat/0611498](https://arxiv.org/abs/cond-mat/0611498), Table I.
A. Kawaguchi, R. Matsui, K. Kobayashi, *The Crystal Structure of Polyethylene at
4.5 K*, Bull. Inst. Chem. Res., Kyoto Univ. **55**(2), 213 (1977), and
A. Kawaguchi, R. Matsui, K. Katayama, *Lattice Distortion of Polyethylene at Low
Temperatures*, ibid. **58**(4) (1980) - both open access in the Kyoto University
Research Information Repository, repository.kulib.kyoto-u.ac.jp.  ACS Omega
**3**, 6244 (2018), [doi:10.1021/acsomega.8b00506](https://doi.org/10.1021/acsomega.8b00506).
J. C. Martin & E. Passaglia, J. Res. NBS **70A**(3), 221 (1966).
C. W. Bunn, Trans. Faraday Soc. **35**, 482 (1939) - the original determination,
never read directly here; every Bunn number above is quoted from a secondary
source that names him.

---

## 2. Chain conformations of the polymorphs

### 2.1 beta = all-trans, alpha = TGTG', gamma = T3GT3G'

**Claim** (`enumerate.py`, `KNOWN_CHAINS["pvdf"]`; DESIGN header and 5.1/5.2
tables): beta is TTTT (all-trans), alpha and delta are TGTG' (`TG+TG-`), gamma
and epsilon are T3GT3G' (`TTTG+TTTG-`).

**Verified.**  Four independent sources agree, and the gamma assignment in
particular is not in dispute:

* Ramer & Stiso: *"The alpha-phase shows a tg+tg- conformation pattern ... The
  beta-phase possesses all-trans conformations ... Regardless of unit cell
  geometry, a tttg+tttg- conformation pattern has been observed [for gamma]
  with eight formula units per unit cell."*
  ([arXiv:cond-mat/0508561](https://arxiv.org/abs/cond-mat/0508561))
* Itoh thesis Table 4.1: Form I `TT`, Form III `T3GT3G-bar`, Form IV `TGT G-bar`,
  Form II `TGT G-bar`.  ([arXiv:1405.5889](https://arxiv.org/abs/1405.5889))
* Su, Strachan & Goddard: *"we indicate the chain conformation with capital
  letters (T means all-T, TG means TGTG' and T3G indicates TTTGTTTG')"*, with
  I = T_p, II = TG_ad, III = T3G_pu, IV = TG_pd.
  ([arXiv:cond-mat/0408156](https://arxiv.org/abs/cond-mat/0408156))
* *"The g-phase of PVDF ... crystallizes in a monoclinic unit cell with four
  chains in a t3g+t3g- conformation"*
  ([arXiv:1706.08068](https://arxiv.org/abs/1706.08068)).

So gamma really is T3GT3G' (= TTTGTTTG'), not something else, and the
`alpha/delta` and `gamma/epsilon` pairings in `KNOWN_CHAINS` are the standard
ones (see 4.2 below).

**Verdict: confirmed.**

### 2.2 beta-PVDF has equal backbone angles and exactly trans torsions

**Claim** (`polymers.py`, the `PVDF` comment; DESIGN 5.8): the accepted Form I
structure has equal C-C-C angles at both backbone carbons and an internal
rotation angle of exactly 180 deg; a DFT study gives 114.4 deg at both carbons,
against 112 deg and a 2.534 A repeat for polyethylene.

**Verified.**  Itoh's Table 4.7 gives, for the Form I crystal, C-C = 1.528 A,
C-F = 1.368 A, C-H = 1.090 A, **C-C-C = 114.4 deg**, F-C-F = 105.4 deg,
H-C-H = 108.0 deg, and the text states *"The calculated internal rotation angle
in Form I crystal is 180 deg"* and *"The lattice constant c and the angle of
C-C-C are slightly, but significantly, larger than those of polyethylene
(2.534 A and 112 deg)."*  There is one C-C-C angle because in the planar zigzag
with the Cm2m coordinates the angles at the CH2 and the CF2 carbon are equal by
the c-translation, which is exactly the point DESIGN 5.8 makes.

Other determinations of the same angle are lower: Ramer & Stiso's own GGA
relaxation gives C1-C2-C1 = 113.6 deg and C1-C2 = 1.53 A, and their Table 2
computes 112.6 deg from Hasegawa's experimental positions and 112.1 deg from
Lando's.  So the range across sources is **112.1-114.4 deg**, with the package's
114.0 deg near the top of it.

**Sources**: Itoh thesis, Tables 4.3 and 4.7 and section 4.3.3
([arXiv:1405.5889](https://arxiv.org/abs/1405.5889)); published as A. Itoh et
al., *Solid-state calculations of poly(vinylidene fluoride) using the hybrid DFT
method: spontaneous polarization of polymorphs*, Polymer Journal **46**, 207-211
(2014), [doi:10.1038/pj.2013.96](https://doi.org/10.1038/pj.2013.96).
Ramer & Stiso, [arXiv:cond-mat/0508561](https://arxiv.org/abs/cond-mat/0508561),
Tables 1 and 2.

**Verdict: confirmed.**

### 2.3 The alternately-deflected beta zigzag is a proposal, not a refinement result

**Claim** (`polymers.py` comment; DESIGN 5.8): Hasegawa's alternately-deflected
zigzag would double the chain repeat to 5.12 A, the intermediate layer line that
doubling requires is absent from the diffraction pattern, and Hasegawa's own fit
used a statistically disordered structure at c = 2.56 A.  The F...F distance
along the chain is 2.56 A against a van der Waals contact of 2.70 A.

**Verified**, essentially sentence for sentence:

> *"In the planar-zigzag structure, the fluorine-fluorine distance is reported
> as 2.56 A, which is equal to the c lattice constant ... This distance is 0.14 A
> less than twice the van der Waals radius of fluorine (2 x 1.35 A = 2.70 A)."*

> *"For their planar-zigzag structure (sigma = 0 deg), Hasegawa et al. determined
> an R of 18.5%.  By varying sigma, they determined a minimum R of 13.5% for
> sigma = 7 deg."*

> *"However, as correctly indicated by Lando et al., these deflections would be
> accompanied by a doubling of the repeat distance along the chain axis
> (2 x 2.56 A = 5.12 A).  Experimentally, an intermediate layer line with that
> period is not present in the diffraction pattern.  Furthermore, the presence of
> these deflections would mean that the assigned space group of Cm2m would no
> longer be correct."*

> *"In response to these contradictions ... Hasegawa et al. propose a
> statistically-deflected zigzag structure ... The statistically-disordered model
> therefore restores the observed periodicity, chain axis length and space
> group."*

The companion paper adds that DFT relaxation collapses the deflection: the
planar zigzag is 0.018 eV per formula unit lower, and the deflection angle
relaxes from the experimental 7 deg to about 3 deg.

**Sources**: N. J. Ramer & K. A. Stiso, *Comparison of polarization and Born
effective charges in alternatively-deflected zigzag and planar-zigzag
beta-poly(vinylidene fluoride)*,
[arXiv:cond-mat/0411564](https://arxiv.org/abs/cond-mat/0411564), section 1;
same authors, [arXiv:cond-mat/0508561](https://arxiv.org/abs/cond-mat/0508561),
section 3.1.

**Verdict: confirmed.**  DESIGN 5.8's withdrawal was correct and is now sourced.

### 2.4 alpha-PVDF's gauche angles are "near +/-45 deg"

**Claim** (DESIGN 2.4): *"real chains deflect (beta-PVDF dihedrals near +/-172
deg, alpha gauche near +/-45 deg)"*.  The beta half of this sentence is already
withdrawn by 5.8 (see 2.2 above: the torsions are exactly trans).

**What the literature says**: the *direction* is supported - Ramer & Stiso note
that the ideal states are t = 180 deg and g+/- = +/-60 deg and that *"the angles
observed experimentally are slightly different"* - but the magnitude is not.
Secondary sources describing Hasegawa's Form II give internal rotation angles of
*179 deg and 45 deg*, but I could not obtain that statement from a source I
could fetch and read, only from search summaries of paywalled text.  Against it,
Itoh's dispersion-corrected DFT of the crystals reports gauche internal rotations
of **59.2 deg (Form III), 61.2 deg (Form IV) and 61.9 deg (Form II)** - that is,
essentially ideal gauche, not 45 deg.

**Verdict: unverified**, with the two candidate values (~45 deg experimental
refinement, ~62 deg DFT) in genuine disagreement.  Do not build an argument on
the 45 deg figure without reading Hasegawa 1972 or Bachmann & Lando,
Macromolecules **14**, 40-46 (1981),
[doi:10.1021/ma50002a006](https://doi.org/10.1021/ma50002a006) directly.
The number does not enter any computation - `polyfind` uses ideal +/-60 deg RIS
states and lets refinement move them - so nothing downstream depends on it.

### 2.5 The real gamma chain's refined torsions

**Claim** (DESIGN 5.2): gamma refines *"to c = 9.27 A (experiment 9.20) with the
deflected trans angles (171-189 deg) and reduced gauche (+/-64 deg) that the
real gamma chain has."*

**Verdict: unverified.**  The refined value of c is `polyfind`'s own output and
the 9.20 A target is confirmed (1.3 above), but I found no published torsion set
for the gamma chain matching "171-189 and +/-64".  The only crystal-structure
torsion figure I could read is Itoh's calculated Form III gauche of 59.2 deg.
The clause "that the real gamma chain has" asserts more than any source I found
supports; it has been softened in DESIGN 5.2.

---

## 3. Spontaneous polarization of beta-PVDF

**Claim** (DESIGN 5.5 table, 5.4 and 5.7): beta-PVDF's polarization is
*"about 0.13 experimentally"*, and the package's 0.1405 C/m^2 lands *"within a
few percent of the experimental value"*.

**What the literature says**: 0.13 C/m^2 is **not** a measurement.  It is the
oldest and crudest *calculated* estimate - the sum of rigid monomer dipoles over
the unit-cell volume.  Itoh's Table 1.1 lists the estimates for Form I:

| Model | P_s (mC/m^2) | Source |
|---|---|---|
| Rigid dipoles | 130 | textbook sum, `P_s = 2 mu_v / abc` |
| Lorentz field | 220 | Kakutani 1970; Mopsik & Broadhurst 1975 |
| Point dipoles | 86 | Purvis & Taylor 1982 |
| Point charges | 127 | Al-Jishi & Taylor 1985 |
| Molecular mechanics | 182 | Carbeck et al. 1995 |
| DFT (Berry phase) | 178 | Nakhmanson et al. 2004 |
| DFT (Berry phase) | 181 | Ramer & Stiso 2005 |
| DFT (Wannier) | 176 | Itoh 2014 |

Itoh's own words: *"As mentioned in Chapter 1, rigid dipole model gives Ps of
ca. 130 mC/m2 for Form I.  The larger value of the calculated Ps suggests that
the polarization is enhanced by 35% in the Form I crystal"*.  Ramer & Stiso:
*"the Ps for the relaxed ... structure of planar-zigzag beta-PVDF was found to be
0.181 C/m2 per unit cell volume.  Our value ... is in excellent agreement with
the value found by a previous DFT study using a similar Berry-phase approach
(0.178 C/m2)"*.

So the modern *calculated* value for a perfect beta crystal is **0.176-0.188
C/m^2**, about 35-45% above the rigid-dipole 0.13.  (A 2021 experimental paper
uses 188 mC/m^2 as "the theoretical limit ... for the neat beta crystal ...
calculated from density functional theory".)

### Measured polarization

What is actually *measured* is the remanent polarization of a poled,
semicrystalline film, and it is smaller and sample-dependent:

| P_r (mC/m^2) | Sample | Source |
|---|---|---|
| 50 | PVDF film, 120 MV/m at 20 C | Furukawa, Date & Fukada, J. Appl. Phys. **51**, 1135 (1980), [doi:10.1063/1.327723](https://doi.org/10.1063/1.327723) - *"the 120-MV/m electric field gave the remanent polarization Pr of 50 mC/m2"* |
| ~65 | stretched PVDF, 240 MV/m | Date, Furukawa & Fukada, J. Appl. Phys. **51**, 3830 (1980), [doi:10.1063/1.328124](https://doi.org/10.1063/1.328124) - *"the remanent polarization is about 65 mC/m2"* |
| 100 | ultradrawn high-crystallinity beta | Nakamura et al., J. Polym. Sci. B **39**, 1371 (2001), [doi:10.1002/polb.1109](https://doi.org/10.1002/polb.1109) |
| 140 (P_s) | biaxially oriented, pure beta | Huang et al., Nat. Commun. **12**, 675 (2021), [doi:10.1038/s41467-020-20662-7](https://doi.org/10.1038/s41467-020-20662-7) |
| 130 +/- 3 | **VDF oligomer**, CF3(CH2CF2)17I, not polymer | Noda et al., J. Appl. Phys. **93**, 2866 (2003), [doi:10.1063/1.1540231](https://doi.org/10.1063/1.1540231) |

Itoh summarises the homopolymer range as *"The reported values of Pr are
scattered in the 50-80 mC/m2 range... The largest value of Pr ever reported is
ca. 100 mC/m2."*

Note the last row: the one place in the literature where **130 mC/m^2 is a
genuine measurement** is a vinylidene-fluoride *oligomer* of extremely high
crystallinity, not PVDF.  That is a plausible route by which the rigid-dipole
0.13 came to be repeated as "the experimental value" - two unrelated 130s.  I
could not find any document that calls 0.13 C/m^2 an experimental value for
beta-PVDF; searches point at one paywalled paper that could not be read, so
that remains unproven either way.

**Verdict: wrong as attributed.**  0.13 C/m^2 is a calculation - the crudest one
in the ladder - not a measurement.  Measured remanent polarizations of PVDF
films are **50-100 mC/m^2** (up to 140 for a biaxially oriented pure-beta film),
and the best calculations for a perfect crystal are 176-188 mC/m^2; a number
near 0.13 sits between the two and belongs to neither.  DESIGN's prose has been
corrected; no computed result changes, because 0.13 only ever appeared as a
comparison target in prose and as one of twelve fit targets in section 5.7.

### 3.1 gamma is polar but weaker than beta; alpha should be zero

**Claim** (DESIGN 5.5 "expectation" column): PE exactly zero; beta about 0.13;
alpha *"should be near zero"*; gamma *"polar, weaker than beta"*.

**Verified** (apart from the 0.13, above).  Itoh computes P_s = 176 mC/m^2 for
Form I, **71 mC/m^2 for Form III (gamma)** and **85 mC/m^2 for Form IV (delta)**,
and notes *"Because samples consisting of only Form III or IV crystal cannot be
obtained for PVDF, Ps for these crystalline forms has not been determined
experimentally"*.  alpha is centrosymmetric and non-polar, so its P_s is zero by
symmetry (4.1 below).

**Verdict: confirmed** for the ordering beta > delta > gamma > alpha = 0.

---

## 4. Which phases are polar and which antipolar

### 4.1 alpha antipolar, beta polar, gamma polar-but-weaker

**Claim** (DESIGN header, 5.2, 5.5).

**Verified**, by three independent sources:

* Ramer & Stiso: *"The chains in alpha-PVDF are packed in an anti-parallel
  manner, yielding a non-polar structure"*; *"The beta-phase possesses all-trans
  conformations ... and parallel chain packing"*; *"a spontaneous polarization in
  beta-PVDF is generated along the b- or y-axis"* (so perpendicular to the chain
  axis c, as DESIGN 5.5 states).
* Itoh: *"molecular packing in crystalline forms I, III and IV is a parallel
  manner; these three polymorphs are polar crystals and are regarded as
  ferroelectrics.  Form II crystal is non-polar, the chains are aligned
  antiparallel."*
* *"Due to the anti-parallel packing of the chains in the unit cell, their dipole
  moments cancel out, rendering the a-phase non-polar and paraelectric.  In
  contrast, the other PVDF polymorphs, i.e. the b-, g- and d-phases, are polar"*
  ([arXiv:1706.08068](https://arxiv.org/abs/1706.08068)).

**Verdict: confirmed.**

### 4.2 delta and epsilon are the polar/antipolar partners of alpha and gamma

**Claim** (DESIGN header; `KNOWN_CHAINS` and `EXPERIMENTAL_CELLS` labels
`alpha/delta` and `gamma/epsilon`).

**Verified.**  *"The d-phase is the polar version of the a-phase; both phases
have the same lattice constants and chain conformation (tg+tg-).  However, in
d-PVDF, every second chain is rotated 180 deg around the chain axis ... the
macromolecules are shifted by half of the c-axis lattice constant"*
([arXiv:1706.08068](https://arxiv.org/abs/1706.08068)).  Ramer & Stiso: *"The
designation of the delta-phase as II_p represents its identification as the polar
form of the non-polar alpha-phase.  It possesses parallel packing of polymer
chains with all its dipoles oriented in the same direction."*  Su et al. use
exactly this pairing (TG_ad = II = alpha, TG_pd = IV = delta, T3G_pu = III =
gamma, T3G_au = the non-polar T3G phase), and record that the non-polar T3G
phase *"was conjectured by Lovinger and proved by Karasawa and Goddard"* - that
is the epsilon phase, the antipolar partner of gamma.

Note that delta shares alpha's cell dimensions but not its space group
(P2_1cn vs P2_1/c, Itoh Table 4.1), which is why the joint `alpha/delta` entry
in `EXPERIMENTAL_CELLS` is legitimate for cell *dimensions* only.

**Verdict: confirmed.**

---

## 5. Relative stability of alpha and beta

**Claim** (DESIGN 5.4): *"the real ordering is nearly degenerate, alpha slightly
favoured."*

**Verified.**  Three sources, agreeing in sign and roughly in magnitude:

* Su, Strachan & Goddard, DFT-GGA, Table I, kcal/mol **per carbon atom**
  relative to the all-trans (beta) crystal: TG_ad (alpha) **-0.59**, TG_pd
  (delta) -0.62, T3G_pu (gamma) -0.51, T3G_au (epsilon) -0.38.  Their text:
  *"our calculations indicate that TG_pd phase is the most stable crystal
  followed by the TG_ad phase."*  Per CH2CF2 monomer that is alpha 1.18 kcal/mol
  below beta.  (Note this is GGA without a dispersion correction.)
  ([arXiv:cond-mat/0408156](https://arxiv.org/abs/cond-mat/0408156))
* Itoh, hybrid DFT (PBE0/cc-pVTZ, dispersion included): *"The relative energy to
  Form I is -2.7 kJ/mol for Form III, -2.5 kJ/mol for Form IV, -2.9 kJ/mol for
  Form II"* per CH2CF2 unit - alpha **0.69 kcal/mol** below beta - followed by
  *"This very small difference of energy may be one of the explanation of the
  stability of all polymorph."*
  ([arXiv:1405.5889](https://arxiv.org/abs/1405.5889))
* *"The a-phase is at ambient conditions the thermodynamically stable
  polymorph"* ([arXiv:1706.08068](https://arxiv.org/abs/1706.08068)).

A fourth source collects the whole spread.  Pelizza & Johnston, *A density
functional theory study of poly(vinylidene difluoride) crystalline phases*,
Polymer **179**, 121585 (2019),
[doi:10.1016/j.polymer.2019.121585](https://doi.org/10.1016/j.polymer.2019.121585)
([open copy](https://strathprints.strath.ac.uk/68678/1/Pelizza_Johnston_Polymer_2019_A_density_functional_theory_study_of_poly_vinylidene_difluoride_crystalline.pdf)),
Table 4, gives *"PVDF phase energies in kJ/mol per monomer relative to the alpha
phase"* - alpha is the reference and positive means less stable.  For beta:
LDA 5.1, PBE 6.5, vdW-DF 3.8, vdW-DF2 2.6, DFT-D2 3.5, with literature rows
PBEsol 3.4, PBE 4.9 (Su et al.), PBE 4.6 (Ranjan et al.), PBE 2.9 (Itoh), and
the MSXX force field at **-3.6** - the one method that gets the sign wrong.
Their earlier paper adds the ordering across the polymorphs: *"all functionals
predict the delta-phase to be the lowest energy polymorph, although the energy
difference of less than 1 kJ/mol is marginal ... The beta-phase has the highest
energy, ranging from 2 kJ/mol for LDA up to 7 kJ/mol for PBE"* (Pelizza, Smith &
Johnston, Eur. Phys. J. Spec. Top. **225**, 1733 (2016),
[doi:10.1140/epjst/e2016-60133-8](https://doi.org/10.1140/epjst/e2016-60133-8)).

Independently of DFT, a free-energy molecular simulation reaches the same
conclusion: *"In concurrence with experimental observations, the alpha-phase was
found to be thermodynamically more stable at normal temperature and pressure
conditions.  The beta-phase was found to be more stable at high and low
temperatures and high pressure"* (Mireja & Khakhar, J. Phys. Chem. B **129**,
6975 (2025), [doi:10.1021/acs.jpcb.5c01058](https://doi.org/10.1021/acs.jpcb.5c01058)).

**Range**: beta lies **2.6 to 6.5 kJ/mol per monomer above alpha**, i.e. alpha
is **0.6-1.6 kcal/mol per monomer** the more stable, with delta essentially
degenerate with alpha (under 1 kJ/mol) and gamma in between.  Pelizza &
Johnston's own summary of why this matters here: *"all phase energies per
monomer vary within a range of less than 10 kJ/mol for all XC functionals.  This
narrow energy range at least partially explains why multiple phases are found
within PVDF samples."*  "Nearly degenerate, alpha slightly favoured" is a fair
summary, and the dispersion correction matters: *"the inclusion of vdW
interactions is, therefore, important to obtain the correct lattice structures,
polarisation and energetics of PVDF polymorphs."*

One thing the sources do **not** say: none of them attributes alpha's prevalence
to kinetics.  alpha is the thermodynamic ground state at ambient conditions in
every method checked, and polymorph coexistence is put down to the smallness of
the gaps.

For contrast, DESIGN 5.4 records that `SimpleFF` puts beta **3.9 kcal/mol per
monomer below alpha** - wrong sign and three to six times too large.  That
remains an accurate statement about the illustrative potential.

**Verdict: confirmed.**

---

## 6. Poly(vinylidene chloride)

### 6.1 PVDC does not adopt a planar zigzag

**Claim** (`docs/CHEMISTRY_EXTENSION.md`, phase 1; repeated in the `PVDC` comment
in `polymers.py`): *"PVDC is known not to adopt the planar zigzag that PVDF's
beta phase does."*

**Verified, and the actual conformation is now on record.**  The crystal
structure is a **glide TGTG' form**, not a planar zigzag:

> *"The molecular conformation is a glide (TGTG') form with internal rotation
> angles of two successive backbone bonds of 175 deg (T) and 49 deg (G'), and
> C-CH2-C and C-CCl2-C bond angles of 123 deg and 114 deg, respectively.  ... The
> deviation of the internal rotation angle from 60 deg (G) to 49 deg (G') and the
> large C-CH2-C bond angle are ascribable to the steric hindrance between the
> chlorine atoms of adjacent monomeric units."*

The cell is monoclinic P2_1 or P2_1/m, a = 6.71 A, b (fibre axis) = 4.68 A,
c = 12.51 A, beta = 123 deg, with four monomeric units (two chains).

**Source**: T. Takahagi, Y. Chatani, T. Kusumoto, H. Tadokoro, *Molecular and
Crystal Structure of Poly(vinylidene chloride)*, Polymer Journal **20**, 883-893
(1988), [doi:10.1295/polymj.20.883](https://doi.org/10.1295/polymj.20.883).

**Verdict: confirmed** - and note the structure independently corroborates the
diagnosis phase 1 offered.  The real chain relieves the Cl...Cl crowding by
opening one backbone angle to 123 deg, which is exactly the degree of freedom
`polyfind` freezes; so the 313 kcal/mol of all-trans strain really is a rigid-
angle artifact, for the reason the document guessed.

### 6.2 PVDC's backbone angles in `polymers.py` are both 114 deg

**Claim** (`polymers.py`, `PVDC`): both backbone atoms get `backbone_angle =
114.0`, with the comment *"Backbone angles are kept equal for the reason given
above"* - the reason being beta-PVDF's equal angles.

**Verified against 6.1**: PVDF's justification does not transfer.  PVDC's angles
are measured as **123 deg at the CH2 carbon and 114 deg at the CCl2 carbon** -
unequal by 9 deg, and unequal for a physical reason.

**Verdict: wrong, and now corrected in the code.**  `polymers.py` carries 123 deg
at the CH2 carbon and 114 at the CCl2 carbon.  The ten-bond all-trans
Lennard-Jones strain falls from 313 kcal/mol to 74, and the energy available by
rotating away from all-trans from 79 kcal/mol to 4.0, so the reference-state
artifact the audit diagnosed is largely gone.  The correction also reproduces the
published *chain*: the torsions that close a TG+TG- repeat under these angles come
out at 175.3 and 49.4 deg with a repeat of 4.677 A, against Takahagi's 175 deg,
49 deg and 4.68 A - and under the old equal angles the same search could only
answer "trans exactly 180" and 4.491 A.  It costs the ability to pack PVDC at
ideal RIS angles, because unequal angles give every repeat 9 deg of curl.  Full
before-and-after in [Changes that were re-measured](#changes-that-were-re-measured).

---

## 7. Other claims

### 7.1 Head-to-head / tail-to-tail defects

**Claim** (DESIGN section 6): *"Head-to-head / tail-to-tail defects (a few
percent in PVDF, known to stabilise beta)."*

**The abundance is confirmed; the effect is disputed.**  A recent measurement
reports defect fractions of *"4.2% at 45 deg C to 5.2% at 90 deg C"* and states
that *"Defect contents between 5% and 10% were reported"* in the literature - so
"a few percent" is right, though the upper end of the quoted range is nearer
10% than "a few".

On the effect, the same paper says the opposite of what DESIGN asserts: *"to
obtain high fractions of beta-phase material the number of defect structures
should be low."*  Against that, the widely repeated result that VDF/TFE
copolymers - which are, structurally, PVDF with extra head-to-head/tail-to-tail
defects - crystallize in the beta form above a few mol% comonomer points the
other way.  I did not find a source I could read that settles it.

**Source**: *Temperature Dependence of the Number of Defect-Structures in
Poly(vinylidene fluoride)*,
[PMC11013231](https://pmc.ncbi.nlm.nih.gov/articles/PMC11013231/).

**Verdict: abundance confirmed (4-6% typical, 5-10% quoted), "known to stabilise
beta" unverified.**  DESIGN's bullet has been rewritten to say so.  Nothing
computed depends on it - the code has no defect placement.

### 7.2 UFF nonbonded parameters

**Claim** (`polymers.py`, `UFF_LJ`): C 3.851/0.105, H 2.886/0.044, N 3.660/0.069,
O 3.500/0.060, F 3.364/0.050, Cl 3.947/0.227, cited to *"UFF: Rappe et al., JACS
1992"*.

**Verified**: all six pairs match the published UFF x1/D1 values for atom types
C_3, H_, N_3, O_3, F_ and Cl exactly, and the citation is right - A. K. Rappe,
C. J. Casewit, K. S. Colwell, W. A. Goddard III, W. M. Skiff, *UFF, a full
periodic table force field for molecular mechanics and molecular dynamics
simulations*, J. Am. Chem. Soc. **114**, 10024-10035 (1992),
[doi:10.1021/ja00051a040](https://doi.org/10.1021/ja00051a040).  Checked against
the parameter table shipped with Open Babel
([UFF.prm](https://raw.githubusercontent.com/openbabel/openbabel/master/data/UFF.prm)),
which cites the same paper.

The derived remark in the `PVDC` comment, *"Chlorine is bulkier than fluorine
(UFF x_i 3.95 vs 3.36 A)"*, follows from those numbers.

**Verdict: confirmed.**

### 7.3 PVDF is "about half crystalline", and the melt is at ~450 K

**Claims** (DESIGN, first sentence; `pipeline.py`, `PipelineConfig.temperature =
450.0  # PVDF melt ~ 450 K`).

**Crystallinity: confirmed.**  *"Since PVDF is a semicrystalline polymer with a
crystallinity of about 50% ..."*
([arXiv:1706.08068](https://arxiv.org/abs/1706.08068)).  Measured crystallinities
in the recent literature cluster in the 40-55% range depending on processing.

**Melt temperature: consistent, but only from secondary sources.**  PVDF's
melting point is widely given as about 170 C (443 K), with individual samples up
to ~178 C (451 K), so an amorphous-ensemble temperature of 450 K is a melt
temperature rather than a sub-Tm one.  I did not find a primary determination
I could read; the number is a modelling choice rather than a claim about a
measurement, and no reported result is sensitive to a few tens of kelvin here.

**Verdict: crystallinity confirmed; 450 K plausible and internally consistent,
not independently sourced.**

### 7.4 The two-chain cell with chain 2 at (1/2, 1/2)

**Claim** (DESIGN 2.3): *"a two-chain cell (PE, alpha-, beta-, gamma-PVDF all
have chain 2 at (1/2, 1/2))"*.

**Mostly confirmed, with one caveat on gamma.**  beta-PVDF is C-centred (Cm2m,
Z = 2 monomers = 2 chains, *"one at the origin and the other at (1/2, 1/2, 0)"* -
Ramer & Stiso), and the gamma cell is C-centred too (Cc, Z = 8 monomers = 2
chains of four monomers).  alpha has Z = 4 monomers = 2 chains.  Polyethylene's
Pnam cell has two chains, the second at (1/2, 1/2).

The caveat: gamma's real packing is usually described as *statistically* parallel
and antiparallel, modelled by a hypothetical **four**-chain cell in space group
C2cm; one source describes the monoclinic cell as containing *"four chains"*.  A
strictly two-chain model of gamma is therefore an idealisation, and `polyfind`'s
"gamma packs parallel" result in DESIGN 5.2 and 5.6 should be read with that in
mind.

**Sources**: Ramer & Stiso
[arXiv:cond-mat/0508561](https://arxiv.org/abs/cond-mat/0508561);
[PMC4379396](https://pmc.ncbi.nlm.nih.gov/articles/PMC4379396/) (*"The chains
pack in a statistical parallel-antiparallel manner, which can be modeled by a
hypothetical four-chain unit cell belonging to space group C2cm"*);
[arXiv:1706.08068](https://arxiv.org/abs/1706.08068).

**Verdict: confirmed for PE, alpha and beta; an idealisation for gamma.**

---

## 8. Chain repeats (the DESIGN 5.1 table)

The "experiment" column of DESIGN 5.1 is the chain-axis repeat of each
polymorph, i.e. the c of the corresponding unit cell, plus one comparison to a
different polymer.

| Chain | quoted | verified | verdict |
|---|---|---|---|
| PE all-trans (2/1) | 2.55 A | see 1.4 | see 1.4 |
| PVDF beta TT | 2.56 A | 2.56 A (Hasegawa; 1.1 above) | confirmed |
| PVDF alpha TG+TG- | 4.62 A | 4.62 A (Hasegawa; 1.2 above) | confirmed |
| PVDF gamma TTTG+TTTG- | 9.20 A | 9.18-9.23 A across sources (1.3 above) | confirmed |
| PE TG helix vs iPP | 6.50 A (iPP) | 6.50 A | confirmed |

### 8.1 Isotactic polypropylene at 6.50 A

**Claim** (DESIGN 5.1, last row): the package's TG polyethylene helix gives
c = 6.35 A, against *"6.50 (iPP)"*.

**Verified**: the alpha form of isotactic polypropylene is monoclinic with
*"a = 6.65 A, b = 20.80 A, c = 6.5 A, beta = 99.8 deg"*, c being the chain axis.
Only the chain-axis **c = 6.50 A** is load-bearing here, and it is 6.50 A in
every source consulted; b and beta vary a little between determinations
(b 20.80-20.96 A, beta 99.3-99.8 deg).

Two further details - four 3/1 helical chains per cell in alternating
handedness, space group C2/c or its subgroup Cc, originally due to Natta and
Corradini - I saw only in search summaries of sources I could not fetch.  They
are standard textbook statements and I have no reason to doubt them, but they
are not verified here and nothing in `polyfind` uses them.

Two honest caveats about the *comparison*, which the row itself already signals
with its "(iPP)":

* it is between different polymers.  A hypothetical all-TG polyethylene helix and
  a real isotactic polypropylene helix share a backbone torsion pattern but not
  bond angles (iPP's backbone angles are wider than polyethylene's), so 6.35
  against 6.50 is a plausibility check on the screw decomposition, not a
  validation against a measurement of the same substance.
* polyethylene does not adopt a TG helix.  The row is testing the geometry
  machinery, not predicting a structure.

**Source** (the sentence actually read): J. Jia and D. Raabe, *Crystallinity and
Crystallographic Texture in Isotactic Polypropylene during Deformation and
Heating*, [arXiv:0811.2412](https://arxiv.org/abs/0811.2412), section 2.2, which
quotes the cell for its own pole-figure analysis.  The primary determination is G. Natta and
P. Corradini, Nuovo Cimento Suppl. **15**, 40 (1960), which I did not read.

**Verdict: confirmed** (the number), with the comparison correctly labelled.

---

## Changes that were re-measured

**All four have now been applied**, together, and every table they touch has been
regenerated: `DESIGN.md` 5.1, 5.2, 5.3, 5.5, 5.7, 5.8, 5.9 and 5.10,
`docs/CHEMISTRY_EXTENSION.md` phases 1 and 3, `docs/VALENCE_FIT.md` section 4 and
its acceptance table, `docs/DFT_FIT.md`'s acceptance table, and
`docs/BENCHMARK.md`'s stage timings.  Each document says which potential produced
its numbers.  The pending list below this section is now empty.

| # | Where | Was | Is | What it did |
|---|---|---|---|---|
| 1 | `polymers.py`, `PVDF.bond_length` | 1.54 A | **1.528 A** (Itoh, Form I crystal) | Applied.  Beta's chain repeat 2.583 -> 2.563 A against an experimental 2.56, error +0.90% -> +0.12%.  **But alpha's and gamma's got worse** (-1.40% -> -2.17% and -0.97% -> -1.74%), and every PVDF density rose about 1% in a potential that already over-predicts them.  Net: better on the polymorph the number comes from, worse on the other two. |
| 2 | `polymers.py`, `PVDC` backbone angles | 114.0 deg at both carbons | **123 deg at CH2, 114 deg at CCl2** (Takahagi et al. 1988; section 6.2) | Applied.  Ten-bond all-trans strain 313 -> 74 kcal/mol, and the drop available by rotating away from trans 79 -> 4.0, so the reference-state artifact is largely gone - the audit's diagnosis, confirmed.  It also **recovers the published chain**: the torsions that close a TG+TG- repeat come out 175.3 and 49.4 deg with c = 4.677 A, against a measured 175, 49 and 4.68.  Cost: with unequal angles no PVDC chain closes at *ideal* torsions (9 deg of curl per repeat), so ideal-angle packing of PVDC is no longer possible. |
| 3 | `pipeline.py`, gamma density | 1.94 | **1.93** for the cell as tabulated (section 1.3) | Applied.  Report-only, but it makes the *reported* agreement worse rather than better: the packed gamma density is 1.960, which was -0.1% against the mismatched 1.94 and is +1.5% against the correct 1.93.  The old agreement was an accident of comparing a density from one determination with a cell from another. |
| 4 | `pipeline.py`, gamma cell | orthorhombic (a, b, c) | monoclinic, **beta ~ 93 deg** (section 1.3) | **Not expressible; documented instead.**  See below. |

### Item 4: the monoclinic cell is not expressible, and this is why

The entry above used to say "the search already fits the cell angle freely, so no
predicted quantity changes".  That is wrong, and checking it was the point of
this exercise.

`polyfind.pack` builds the cell as `a = (a, 0, 0)`, `b = (b cos g, b sin g, 0)`
and `c = (0, 0, c)`, with the chain axis along z.  The free angle is therefore the
crystallographic **gamma**, between a and b, in the plane perpendicular to the
chains; `default_bounds` pins it at 90 deg unless `pack(..., gamma_free=True)` is
passed, and with that flag it does move (gamma-PVDF settles at 87.2 deg, 0.08
kcal/mol per monomer below the right-angled cell).  **alpha and beta are 90 deg by
construction and no flag changes that**: nothing in the parametrisation can tilt a
lattice vector out of the plane normal to the chain axis, and `to_cif` writes both
as 90 unconditionally.

gamma-PVDF's unique angle is beta = 93 deg, between a and the chain axis c.  It is
exactly the one the model cannot represent.  The package therefore **approximates
the gamma cell as orthorhombic**, and `pipeline.REFERENCE_CELL_ANGLES` now records
the real angles and says so, rather than the table silently implying right angles.

What the approximation costs, measured: a lattice translation along a would carry
a z-component of `a cos(beta)` = -0.26 A, which the model must set to zero, and
the cell volume changes by `sin(beta)` = 0.9986, so the best density the model
could reproduce differs from the monoclinic one by 0.14%.  That is well inside the
0.5% granularity of a two-decimal density and far inside this potential's own
3 to 9% cell-edge errors, so no reported number is affected.  It would matter to a
structure comparison at crystallographic precision, and to nothing here.

### Not re-measured, and why

Two things these corrections moved that have deliberately *not* been regenerated:

* **the fits themselves.**  `pvdf-crystal-fit`, `pvdf-dft-fit` and
  `pvdf-dft-valence` were fitted at the old geometry, over 25 minutes, 40 seconds
  and 260 seconds respectively.  Their *parameters* are unchanged, so what the
  documents quote as fit outcomes (objective values, held-out errors, ablation
  tables) still describes those runs.  Everything that is a straight *evaluation*
  of a shipped preset - acceptance tests, polarizations, cell predictions - has
  been re-measured and is marked as such.  Refitting is a separate job.
* **the inline records in `src/polyfind/fitting.py`.**  That file is not owned by
  this change and was not edited, so its comments still quote the pre-correction
  acceptance numbers (-4.54 and -2.96 kJ/mol per monomer, RIS margin 0.18, PVDC
  strain 313).  The re-measured values are in `docs/DFT_FIT.md` and
  `docs/VALENCE_FIT.md`; the comments should be brought into line the next time
  that file is touched.

## Changes that need re-measurement

None outstanding.  Two numbers that were candidates for the list above and are
deliberately *not* on it:

* the polyethylene c of 2.55 A - an intermediate draft of this file had it down
  as too high, on the strength of Bunn's 1939 value of 2.534 A.  That was wrong:
  the modern room-temperature determinations give 2.550-2.551 A, so 2.55 is
  correct and there is nothing to change (section 1.4);
* the PVDF backbone angle of 114.0 deg - within the 112.1-114.4 deg spread of
  the determinations, and the value the package uses is the one the DFT study it
  cites reports (section 2.2).
