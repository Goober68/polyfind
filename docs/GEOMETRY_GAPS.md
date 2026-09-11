# Backbone geometry for the non-PVDF chemistries: what exists and what does not

The screen (`docs/SCREEN.md`) explicitly declined to stand behind any result for
CFE or CDFE because both borrow PVDF's equal 114-degree backbone angles. The
precedent for fixing that is strong: correcting PVDC's angles to the measured 123
at CH2 and 114 at CCl2 dropped its all-trans strain from 313 kcal/mol to 74 and
recovered the published glide chain's torsions to within half a degree. That was a
literature lookup.

This is the same lookup for the rest. **Nothing here is applied**; it is queued the
way the four earlier corrections were, so the documented tables are regenerated
once. Two of the four entries are not parameter gaps at all but limitations worth
knowing.

## AN / polyacrylonitrile: there is no reference geometry to borrow

This is the important entry. Polyacrylonitrile does not have a settled
crystallographic repeat. The literature describes it as pseudo-crystalline: a
helical conformation has been proposed with four monomers per repeat and an
orthorhombic cell, a syndiotactic configuration has been invoked to explain the
meridional reflections, and molecular-dynamics work finds that **most
configurations, regardless of tacticity, lose periodicity along the chain axis on
relaxation**, giving pseudo-crystalline rather than crystalline structures.

That is a limitation of our method applied to AN, not a missing number. This
package's entire conformational search assumes a crystallographic repeat exists
and enumerates periodic sequences; if the real polymer does not maintain one, then
the AN column of the screen is answering a question the material does not pose. It
should carry that caveat regardless of what angles we use.

It also connects to something Sarco measured independently: their accepted
VDCN/VDF copolymer endpoint has sixteen of 192 dihedrals more than thirty degrees
from any rotational-isomeric state, about one kink per chain. Loss of chain
periodicity in nitrile-bearing backbones appears from two directions at once.

## CFE and CDFE: no homopolymer structure; one related datum

Neither exists as a characterised homopolymer crystal - both are copolymer
comonomers in practice - so there is nothing to correct *to* directly. The closest
characterised analogue is poly(chlorotrifluoroethylene), a fully halogenated
backbone:

| quantity | PCTFE | PVDF (for comparison) |
|---|---|---|
| C-C-C backbone angle | **118 deg** | 114.4 deg |
| C-C bond | 1.57 A | 1.528 A |
| C-F bond | 1.34 A | 1.35 A (ours) |
| C-Cl bond | 1.74 A | 1.77 A (our PVDC) |
| F-C-F | 109.47 deg | 105.4 deg |
| chain | 17/1 helix, hexagonal a = 6.37 A | planar zigzag |

So a fully halogenated backbone opens to 118 degrees against PVDF's 114.4, and its
C-C bond is longer. That is a direction rather than a value for CFE and CDFE, and
the PVDC result adds a pattern worth noting: there, the angle that opened was the
**CH2** one, to 123 degrees, while the substituted carbon stayed near 114. The
crowded partner pushes the unsubstituted carbon open rather than widening itself.

Applying that pattern to CFE (CH2 with CFCl) or CDFE (CHCl with CF2) would be
interpolation between PVDF, PVDC and PCTFE, not a measurement, and should be
labelled as such if anyone does it. It would still likely beat the current
borrowed 114/114, since that is the one combination the three references agree is
wrong for a halogen-bearing backbone.

## VDCN and FANOME

Not searched. Poly(vinylidene cyanide) is not a stable isolable homopolymer, so a
crystal structure is unlikely to exist, and FANOME is already registered but
deliberately unfitted because its all-trans geometry is unphysical in this model.

## Recommendation

Three things, in order of value.

1. **Add the AN caveat now**, in the screen and in the polymer definition. It costs
   nothing and it is the difference between a number a reader trusts and one they
   do not. A polymer whose chains do not maintain periodicity is outside what this
   method models, and saying so is more useful than any angle we could pick.
2. **Leave CFE and CDFE at 114/114 until someone decides** whether an interpolated
   geometry labelled as interpolation is better than a borrowed one labelled as
   borrowed. Our view is that it is, given three independent references agree the
   current value is wrong in a known direction, but it is a judgement rather than a
   correction and should be made deliberately.
3. If the interpolation is wanted, batch it with the other queued changes so the
   tables regenerate once, following `docs/REFERENCES.md`'s existing list.

Sources: [PAN crystal structure](https://www.sciencedirect.com/science/article/abs/pii/0014305774901475),
[PAN crystalline and pseudo-crystalline phases from molecular dynamics](https://www.sciencedirect.com/science/article/abs/pii/S0032386118308620),
[PAN hierarchical structure](https://pmc.ncbi.nlm.nih.gov/articles/PMC9064467/),
[PCTFE crystallization](https://www.nature.com/articles/pj197216),
[PTFE and PCTFE structure refinement](https://arxiv.org/pdf/1309.6943).
