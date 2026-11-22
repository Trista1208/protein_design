# What each stage does — in plain language

Audience: a teammate who knows biology but not this toolchain, or a judge who
wants to know what we did without reading code.

The full technical reasoning lives in `../HANDOVER_kemp_eliminase.md`. This file
explains *what* and *why*; the user guide explains *how to run it*.

---

## The task in one paragraph

We must invent a brand-new protein that speeds up one specific chemical reaction
(a **Kemp elimination**). We are given a tiny fragment of a known enzyme's active
site — three pieces: an **aspartate** that does the chemistry, a **glutamine**
that steadies the molecule mid-reaction, and the **target molecule itself**. Our
job is to build a whole protein around those three pieces. The catch: it must not
resemble any published enzyme, or it is disqualified.

Think of it as being handed three gears in a precise arrangement and being told
to design an entirely original machine that holds those gears exactly so.

---

## The two chemicals, and why people confuse them

| Molecule | Role | Where it appears |
|---|---|---|
| **5-nitrobenzisoxazole** | the real **substrate** — what the enzyme actually consumes | the reaction we're accelerating |
| **6-nitrobenzotriazole** (`6NT`) | a **stand-in** that looks like the molecule mid-reaction | what we score against |

We score against the stand-in because it mimics the fleeting half-way state that
an enzyme must grip tightly. Two near-identical decoys exist and the code refuses
them by name: `H5J` (the real substrate, from a different structure) and
5-nitrobenzotriazole (nitro group in the wrong position).

---

## Why the supplied files contain a trap

The organisers cropped our fragment out of a real crystal structure, **6TU6**.
That structure is a *deliberately broken* version of the enzyme — the glutamine
was mutated away to measure how much it contributes. Someone then **computationally
put the glutamine back** before handing it to us.

How we know: the aspartate and the target molecule have realistic, refined
measurements. The glutamine has a B-factor of exactly **0.00 on all 17 atoms**
and no anisotropic record — the signature of an atom that was typed in, not
observed.

**Consequence:** the glutamine's position is somebody's *prediction*, not a
measurement. Our scoring treats the aspartate contact as trustworthy (60% of the
geometry weight) and the glutamine contact as a weak hint (15%, explicitly
flagged in every report).

### The second trap: two files, two coordinate systems

We get two files. Their atoms sit in **different places in space** while
describing the *same* arrangement — like the same photo printed at two positions
on a page. Our audit proves this: after rotating one onto the other they agree to
**0.0005 Å**, but as supplied they are **68.7 Å apart**.

One file carries an extra marker (`ORI`) telling the design software where to
centre the new protein. Copy that marker's numbers into the *other* file's
coordinate system and the origin lands ~70 Å from the active site. **The software
does not complain.** It quietly builds proteins centred on empty space, and every
design in that batch is rubbish.

The code refuses to guess which coordinate system you meant, and asserts the
marker sits within 12 Å of the target molecule before every single run.

---

## The free experiment the organisers handed us

This is the most valuable idea in our plan.

Because the fragment came from the *broken* enzyme, both versions are published:

- **HG3.17** — the working enzyme, glutamine present. Very good at the reaction.
- **HG3.17-Q50M** — the same protein with the glutamine deleted. Measurably worse.

So before designing anything, we run **both** through our complete scoring
pipeline. If our scoring cannot tell that the working one is better, then our
scoring is not measuring the thing that matters, and **every number we produce
afterwards is noise**.

This costs about an hour. The pipeline **refuses to run the large design stage**
until this check passes — that refusal is enforced in code, not a note in a
README.

It also gives us a **noise floor**: how much our own scores wobble between runs
on a protein we *know* is real. Any two designs closer together than that wobble
are indistinguishable, and ranking within it is superstition.

---

## The milestones

### M0 — Audit and controls *(CPU, ~2 h)*
Read the supplied files safely, prove the two coordinate systems story, verify
our measurements match the reference table, strip the marker atom, and write one
clean file everything else uses. Then score the two control enzymes.

**Why it matters:** three of the four ways this project could silently fail are
caught here.

### M1 — Design smoke test *(GPU, ~2 h)*
Run the backbone generator once, on the smallest possible setting, purely to
measure how much graphics memory and wall-clock time it needs. The documentation
does not state its memory requirement, so we measure it.

### M2 — Sequence design and folding smoke test *(GPU, ~2 h)*
Given a backbone shape, decide which amino acids to put at each position
(**LigandMPNN**), then independently predict what that sequence would actually
fold into (**Boltz-2**) to check it matches the intended shape.

**Critical detail:** the sequence designer will happily mutate our catalytic
aspartate into something else, destroying the enzyme. We pin those positions and
then *assert* they survived, failing loudly if not.

### M3 — Scoring and eligibility *(CPU, ~3 h)*
Four official measures:

| Metric | Plain meaning |
|---|---|
| **iPDE** | how confident the folding model is about the protein–molecule interface |
| **Geometry** | do the catalytic pieces end up in the right arrangement? |
| **Base burial** | is the catalytic aspartate properly enclosed in a pocket, not on the surface? |
| **Consistency** | do repeated runs agree, or is one lucky result? |

Plus three hard disqualifiers: over **260 amino acids**, over **70% sequence
identity** to a published enzyme, or a structural match (**TM-score > 0.7**).

Two honest gaps, both deliberate rather than papered over:
- The exact field behind "iPDE" is unconfirmed, so we **read** it from the
  folding output and never invent a formula. Absent field → reject the candidate.
- The "within 12 Å" burial measurement could be anchored on three different
  atoms. We compute **all three** and report them side by side.

### M4 — Scale *(GPU, overnight)*
Generate many backbones across several active-site arrangements and lengths,
design sequences for each, filter hard, and run the expensive 5-replicate folding
**only** on survivors. Blocked until M0's controls pass.

### M5 — Select the final five *(CPU, ~2 h)*
Pick five designs that are **different from each other**, not five variations on
one pocket. Evidence from a comparable competition: of 400 proteins tested, 95%
could be produced but only 14% worked. Diversity beats depth.

Two tie-breakers:
- **ESM2 likelihood** — does the sequence look like something a real organism
  might make? Proxy for "will actually express in the lab". The official scoring
  rewards nothing about this, and five designs that don't express are five empty
  wells.
- **FieldTools** — measures the electric field the protein projects onto the
  breaking bond. A genuinely different kind of evidence from folding confidence.
  Reported as a **proxy**, never as an activation energy.

---

## Five different active sites, not one

Designing five variations on one pocket is the most likely way to waste all five
laboratory slots. We build at least three of:

1. **Canonical** — aspartate + glutamine, exactly as supplied.
2. **Crystallographic glutamine** — the real observed rotamer rather than the modelled one.
3. **Backbone amide donor** — no glutamine at all; the protein's own backbone does the steadying.
4. **Alternative donor** — serine/threonine/tyrosine instead of glutamine.
5. **Different base** — glutamate instead of aspartate (reaches further).

---

## What we deliberately did NOT do

| Rejected | Why |
|---|---|
| Improving an existing Kemp eliminase | Disqualified by the novelty rules by construction. Kept only as control sequences. |
| Full AI.zymes evolutionary loop | Needs a compliant starting scaffold we don't have until M4, plus heavy setup. We kept its **FieldTools** component, which is the genuinely novel idea. |
| denovoKemp / Rosetta Matcher | A multi-month CPU campaign in the original paper. Not reproducible on a hackathon clock. We kept two insights: aim for ~190–250 residues, and don't assume the canonical arrangement is optimal. |
| Web services (FuncLib, PROSS) | Queued external services with unpredictable turnaround. |

The winning entry in the comparable competition won by optimising an existing
protein — **a route these rules deliberately close.** Recognising that early is
why effort went into one well-executed de novo path instead of four half-finished ones.
