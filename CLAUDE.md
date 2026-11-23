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
- Both are read against the ceiling the heat budget sets whatever the chemistry: force density
  times speed is at most the heat budget over pi times the mechanical loss tangent. With the
  tables' settings (4 MW/m3, tan delta mechanical 0.05) that is 25 W/cm3. Stiffness and drive
  frequency trade force for speed along it; only the heat budget and the mechanical loss move it.
- The modulus (0.2 GPa), the mechanical loss tangent (0.05) and the heat budget (4 MW/m3) are
  assumed, the same for every candidate, and not measured or computed here. Say so wherever a
  force density or a speed is given, until the PI supplies measured ones.

