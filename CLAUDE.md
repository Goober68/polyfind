# Project memory

Standing instructions from the PI that hold across sessions. The state of the work is in
`RESUME.md`; this file is for what must not be forgotten when that state changes.

## The objective is a point on the strain / tan delta / frequency curve (PI, 2026-11-14)

The PI's words: "strain is NOT the sole metric by which to gauge. What we want is the to be able
to evaluate candidates for the optimal point on the strain/tanD/freq curve. For example: a high
frequency, low loss moderate strain polymer will outperform our CFE."

Scope: every evaluation, ranking and recommendation of a candidate material in this project.

How it is applied (an interpretation, not the PI's words):

- A candidate is reported on all three axes together: strain at the operating point, tan delta,
  and the frequency it follows. It is ranked on the combination, never on strain alone.
- A candidate with less strain than the CFE terpolymer is not out. If it is faster or less lossy
  it can be the better material, and the write-up has to say so.
- A result that reads only on strain (a polar fraction at the bias, a margin's place in the
  switching window) is a partial result until its loss and its rate are beside it.

## Blocked force is the fourth part (PI, 2026-11-15)

The PI's words: "Do we have a way to estimate blocked force?  Thats the 4th part of the equation"

Scope: the same as above. Every evaluation, ranking and recommendation of a candidate material.

How it is applied (an interpretation, not the PI's words):

- A candidate is reported on four things together: strain, tan delta, the frequency it follows,
  and the blocked stress it can hold (force per area; a force needs the actuator's cross-section).
- Blocked stress is stiffness times free strain, so it is a separate axis only where the
  candidate's own stiffness is known. With one assumed modulus for every candidate it is the
  strain axis again, and the write-up has to say which of the two it is.

The PI's words on its unit (2026-11-15): "Actuator has no fixed xsection, so establish an
abstract, say force = N/nm2"

How it is applied (an interpretation, not the PI's words):

- Blocked force is reported per area, never as a force for an assumed cross-section.
- The unit is pN/nm2, which is the same number as MPa and N/mm2. Plain N/nm2 would put the
  values near 1e-13; say so once where the unit is introduced.

## Speed is the fifth part, per length (PI, 2026-11-23)

The PI's words: "we need a mm/s or cm/s number to go along with our pN/nm2 i guess"

Scope: the same as above. Every evaluation, ranking and recommendation of a candidate material.

How it is applied (an interpretation, not the PI's words; the PI was asked which normalisation
and which definition and has not said, so these are the defaults until told otherwise):

- Speed is reported per unit of active length, mm/s per mm (the same number in cm/s per cm, and
  numerically 1/s): a speed needs the actuator's length as a force needs its cross-section.
- It is the loss model's stroke rate, AC strain per cycle times drive frequency
  (`polyfind.loss`: `stroke_rate_per_lambda`). A peak sinusoidal velocity is 2 pi times it; say
  which one a number is.
- Force per area times speed per length is power per volume: pN/nm2 times mm/s per mm is W/cm3.

## The aim is a step up in force density and speed (PI, 2026-11-23)

The PI's words: "Looking for a step up in force density and speed". Said after the report of
`docs/MLIP_SECOND_OPINION.md` 3.26, and after: "Im shocked and troubled that 6.6% better 8s all
weve come up with - its like all were looking for is the exact same activation shape in different
materials".

Scope: what a candidate is for. It sits above the lists of axes: those say what is reported, this
says what counts as a result.

How it is applied (an interpretation, not the PI's words):

- A candidate that moves frequency or tan delta and leaves force density and speed where the
  film's composition has them is a partial result, and the write-up says so in those words.
- Force density and speed are reported on ONE modulus and under ONE heat budget. A strain that the
  heat budget allows a soft film is not multiplied by a crystal's stiffness: a film that stiff
  would not be allowed that strain. (The "blocked force through the crystal" of 3.20 to 3.26 does
  this and overstates; see RESUME.)
- (WITHDRAWN 2026-11-23, the same day it was written: a "ceiling" of heat budget over pi times
  the mechanical loss tangent on force density times speed. It was this file's interpretation,
  not the PI's, and the PI's next two messages showed what it rests on: see the next section.)
- The modulus (0.2 GPa), the mechanical loss tangent (0.05) and the heat budget (4 MW/m3) are
  assumed, the same for every candidate, and not measured or computed here. Say so wherever a
  force density or a speed is given, until the PI supplies measured ones.

## No candidate is pinned to 34 kHz (PI, 2026-11-23)

The PI's words, on a ranking that put every candidate at 34 kHz: "Why the fuck do i care about
34khz - what is the operating frequency of these new ones"

Scope: every table, ranking and speed number of a candidate.

How it is applied (an interpretation, not the PI's words; the PI was asked what frequency the
device needs and has not said):

- 34 and 66 kHz are the band the programme's own optimiser assumed, carried into
  `examples/pin_design.py` on 2026-11-04. They are not a requirement the PI gave here, and a
  candidate is not ranked at them as if they were.
- A candidate is reported at the frequency it can run at, with what sets that frequency: a
  switch's relaxation, or, for a lattice that only answers linearly, the part's own mechanical
  resonance and its heat. Say which, and say when the limit is not computed.
- Speed is given so that it can be read at any frequency (per MHz, or at the candidate's own
  limit), not at one borrowed drive.

## Free motion stores and returns energy; work and heat begin at the blocking force (PI, 2026-11-23)

The PI's words: "the energy dissipated depends on the work done not your theoretical numbers - the
heat bidget you have assumes some efficiency level?  What is it?" and then: "The reason i make the
point about work done is because non blocked motion is like charging a battery, the potential is
stored and released with a cycle efficiency of maybe 80%.  When blocking force is encountered, then
energy starts being converted into work/heat."

Scope: how dissipation, efficiency, heat and work are reckoned for any candidate, and so every
force density, speed and power number.

How it is applied (an interpretation, not the PI's words):

- Two regimes, kept apart. Free (unblocked) motion: energy is stored and returned each cycle, and
  what is lost is the round trip's shortfall (the PI's figure: about 80% returned). Motion against
  a force: stored energy is converted to work and heat.
- `polyfind.loss.operating_point` does NOT do this. It charges the mechanical loss at the full
  free strain whether or not work is being done, and credits `Y S^2 / 4` of work a cycle as if
  always at a matched load. An efficiency, a heat-limited strain or a "ceiling" taken from it
  pairs loaded work with unloaded loss. Until the model is rebuilt on the two regimes, say this
  beside any such number and do not present it as a limit of the material.
- The mechanical loss tangent (0.05), the film modulus (0.2 GPa) and the heat budget (4 MW/m3) are
  assumed. The PI's 80% round trip is the only figure of this kind that comes from the device; ask
  what it was measured on before using it as a number.

