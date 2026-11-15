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
