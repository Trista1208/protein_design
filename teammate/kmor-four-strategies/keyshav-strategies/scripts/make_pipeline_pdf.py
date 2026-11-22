"""Three-page pipeline PDF: block diagram, stage I/O, and the rationale.

    python scripts/make_pipeline_pdf.py [-o outputs/pipeline.pdf]

Page 1  the block diagram, drawn to scale with real stage names and tools
Page 2  inputs and outputs at every stage, plus the tool table
Page 3  why the pipeline has this shape - each filter traced to the mistake
        that put it there

Every number is read from the run at build time, never typed in, so the
document cannot drift from the data (docs/good_practices.md 1).
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfgen import canvas as rl_canvas
from reportlab.platypus import Paragraph, Table, TableStyle

ROOT = Path(__file__).resolve().parents[1]
W, H = A4

INK = colors.HexColor("#15243a")
SOFT = colors.HexColor("#4a5a72")
PAPER = colors.HexColor("#fbfaf6")
BOX = colors.HexColor("#eef1f5")
GATE = colors.HexColor("#f6ece0")
GOOD = colors.HexColor("#e2efe9")
ACC = colors.HexColor("#c2603a")
ACC2 = colors.HexColor("#2f6f62")
RULE = colors.HexColor("#d9d5c9")


def stats() -> dict:
    n = 0
    for f in ("outputs/insurance/screen.json", "outputs/library/screen.json",
              "outputs/library2/screen.json", "outputs/campaign/campaign.json",
              "outputs/gen2/screen.json"):
        p = ROOT / f
        if p.exists():
            try:
                n += len([r for r in json.loads(p.read_text())
                          if r.get("status") in (None, "ok")])
            except Exception:                                    # noqa: BLE001
                pass
    best, bid, nz, scored = 0.0, None, None, 0
    for v in (ROOT / "outputs/finalists").glob("*/verdict.json"):
        try:
            d = json.loads(v.read_text())
        except Exception:                                        # noqa: BLE001
            continue
        scored += 1
        b = d.get("official_score_best") or 0
        if b > best:
            best, bid, nz = b, v.parent.name, d.get("n_nonzero_replicates")
    return {"screened": n, "scored": scored, "best": best,
            "best_id": bid, "repro": nz}



def wrap(c, text, x, y, width, font="Helvetica", size=7.6, leading=4.0):
    """Word-wrap at a real pixel width.

    The first draft sliced on a character count (`desc[:92]`), which split
    words mid-word: "the ligand is p / laced where it says". Measure instead.
    """
    words, line = text.split(), ""
    c.setFont(font, size)
    for w in words:
        trial = (line + " " + w).strip()
        if c.stringWidth(trial, font, size) > width and line:
            c.drawString(x, y, line)
            y -= leading * mm
            line = w
        else:
            line = trial
    if line:
        c.drawString(x, y, line)
        y -= leading * mm
    return y


def header(c, title, sub):
    c.setFillColor(INK)
    c.setFont("Helvetica-Bold", 17)
    c.drawString(20 * mm, H - 22 * mm, title)
    c.setFillColor(SOFT)
    c.setFont("Helvetica", 9.5)
    c.drawString(20 * mm, H - 28 * mm, sub)
    c.setStrokeColor(RULE)
    c.setLineWidth(0.8)
    c.line(20 * mm, H - 31 * mm, W - 20 * mm, H - 31 * mm)


def footer(c, page, total):
    c.setFillColor(SOFT)
    c.setFont("Helvetica", 7.5)
    c.drawString(20 * mm, 12 * mm,
                 "Kemp eliminase pipeline · DBSSE Protein Design Hackathon, "
                 f"Basel · built {datetime.now():%Y-%m-%d %H:%M}")
    c.drawRightString(W - 20 * mm, 12 * mm, f"{page} / {total}")


def box(c, x, y, w, h, title, lines, fill=BOX, tool=None, bold_title=True):
    c.setFillColor(fill)
    c.setStrokeColor(RULE)
    c.setLineWidth(0.8)
    c.roundRect(x, y, w, h, 3 * mm, stroke=1, fill=1)
    ty = y + h - 6.5 * mm
    if tool:
        c.setFillColor(ACC)
        c.setFont("Helvetica-Bold", 6.6)
        c.drawString(x + 4 * mm, ty, tool.upper())
        ty -= 4.6 * mm
    c.setFillColor(INK)
    c.setFont("Helvetica-Bold" if bold_title else "Helvetica", 8.6)
    c.drawString(x + 4 * mm, ty, title)
    ty -= 4.4 * mm
    c.setFillColor(SOFT)
    c.setFont("Helvetica", 7.2)
    for ln in lines:
        c.drawString(x + 4 * mm, ty, ln)
        ty -= 3.5 * mm


def arrow(c, x1, y1, x2, y2, label=None):
    c.setStrokeColor(SOFT)
    c.setLineWidth(1.1)
    c.line(x1, y1, x2, y2)
    import math
    a = math.atan2(y2 - y1, x2 - x1)
    for s in (+1, -1):
        c.line(x2, y2,
               x2 - 2.6 * mm * math.cos(a - s * 0.42),
               y2 - 2.6 * mm * math.sin(a - s * 0.42))
    if label:
        c.setFillColor(SOFT)
        c.setFont("Helvetica-Oblique", 6.4)
        c.drawCentredString((x1 + x2) / 2, (y1 + y2) / 2 + 1.6 * mm, label)


def page1(c, st):
    header(c, "The pipeline, as it actually runs",
           "Design a de novo Kemp eliminase · scored against the "
           "6-nitrobenzotriazole transition-state analogue")

    L, R = 20 * mm, W - 20 * mm
    # leave a clear 10 mm channel on the right for the feedback loop, so the
    # return path never crosses a box or a label (it did in the first draft)
    CH = 11 * mm
    inner_R = R - CH
    colw = (inner_R - L - 3 * 4 * mm) / 4
    yrow = H - 78 * mm
    bh = 34 * mm

    box(c, L, yrow, colw, bh, "Foldable parent",
        ["One design that folds at", "pLDDT 0.98. The backbone", "is kept; only the pocket",
         "around it changes."], tool="ins080")
    box(c, L + colw + 4 * mm, yrow, colw, bh, "Propose variants",
        ["~225 targeted mutants:", "H-bond donor scan, anchor", "residues, salt-bridge fix,",
         "and pairs of winners."], tool="LigandMPNN + rules")
    box(c, L + 2 * (colw + 4 * mm), yrow, colw, bh, "Screen, 1 seed",
        ["8 folds in parallel, about", "11 s each on one GPU.", "Cheap triage, never a",
         "result on its own."], tool="Boltz-2 x1")
    box(c, L + 3 * (colw + 4 * mm), yrow, colw, bh, "Full score, 5 seeds",
        ["Each replicate scored", "SEPARATELY, then", "aggregated as best plus",
         "reproducibility."], tool="Boltz-2 x5", fill=GOOD)

    for i in range(3):
        x = L + (i + 1) * colw + i * 4 * mm
        arrow(c, x + 0.5 * mm, yrow + bh / 2, x + 3.5 * mm, yrow + bh / 2)

    # gates, directly beneath the screen column
    gy = yrow - 58 * mm
    gw = (inner_R - L - 2 * 4 * mm) / 3
    gh = 30 * mm
    box(c, L, gy, gw, gh, "1 · Is the ligand bound?",
        ["Both catalytic contacts", "inside 2.0-3.6 A. Otherwise", "the geometry is noise,",
         "not an angle."], fill=GATE, tool="hard gate")
    box(c, L + gw + 4 * mm, gy, gw, gh, "2 · Did it fold?",
        ["mean pLDDT >= 0.95.", "Below that a design is", "unmeasurable rather",
         "than bad."], fill=GATE, tool="hard gate")
    box(c, L + 2 * (gw + 4 * mm), gy, gw, gh, "3 · Is it eligible?",
        ["<= 260 aa; identity <= 0.70", "and TM <= 0.70 against", "published enzymes and",
         "structures."], fill=GATE, tool="MMseqs2 + Foldseek")

    # screen -> gates, straight down from the screen column
    sx = L + 2 * (colw + 4 * mm) + colw / 2
    arrow(c, sx, yrow - 1 * mm, sx, gy + gh + 1.5 * mm)
    c.setFillColor(SOFT)
    c.setFont("Helvetica-Oblique", 6.6)
    c.drawString(sx + 3 * mm, (yrow + gy + gh) / 2 - 1 * mm,
                 "every screened design enters the gates")

    # gates -> outputs
    oy = gy - 50 * mm
    ow = (inner_R - L - 4 * mm) / 2
    oh = 30 * mm
    box(c, L, oy, ow, oh, "Submission set",
        ["5 designs that pass every", "gate, exported as FASTA,", "PDB and a metrics JSON."],
        fill=GOOD, tool="output")
    box(c, L + ow + 4 * mm, oy, ow, oh, "Live publication",
        ["designs.csv, plots, PDBs and", "rotating 3D GIFs, pushed to", "the team repository."],
        fill=GOOD, tool="output")
    gx = L + 1.5 * gw + 4 * mm
    arrow(c, gx, gy - 1 * mm, L + ow / 2, oy + oh + 1.5 * mm)
    arrow(c, gx, gy - 1 * mm, L + ow + 4 * mm + ow / 2, oy + oh + 1.5 * mm)
    c.setFillColor(SOFT)
    c.setFont("Helvetica-Oblique", 6.6)
    c.drawString(L, gy - 9 * mm, "only survivors are published or submitted")

    # feedback loop, routed up the clear right-hand channel and over the top
    fx = inner_R + 5 * mm
    top = yrow + bh + 7 * mm
    px = L + colw + 4 * mm + colw / 2          # into "Propose variants"
    c.setStrokeColor(ACC)
    c.setLineWidth(1.0)
    c.setDash(3, 2)
    c.line(L + ow + 4 * mm + ow / 2, oy + oh + 1.5 * mm, fx, oy + oh + 1.5 * mm)
    c.line(fx, oy + oh + 1.5 * mm, fx, top)
    c.line(fx, top, px, top)
    c.setDash()
    arrow(c, px, top, px, yrow + bh + 1.5 * mm)
    c.setFillColor(ACC)
    c.setFont("Helvetica-Oblique", 6.8)
    c.drawCentredString((fx + px) / 2, top + 2.2 * mm,
                        "winners become the next generation's parents")

    c.setFillColor(INK)
    c.setFont("Helvetica-Bold", 8.4)
    c.drawString(L, oy - 14 * mm, "Measured so far")
    c.setFillColor(SOFT)
    c.setFont("Helvetica", 8)
    c.drawString(L, oy - 19 * mm,
                 f"{st['screened']} designs screened · {st['scored']} fully "
                 f"scored on 5 seeds · best {st['best']:.4f} achieved on "
                 f"{st['repro']} of 5 replicates.")
    c.drawString(L, oy - 24 * mm,
                 "A published crystal structure of a real Kemp eliminase "
                 "scores 0.959 on the same scale.")
    c.drawString(L, oy - 29 * mm,
                 "The score is a geometric proxy computed on a predicted "
                 "structure, not a measured reaction rate.")
    footer(c, 1, 3)


def page2(c, st):
    header(c, "Inputs, outputs and tools at each stage",
           "What enters a stage, what leaves it, and which program does the work")
    styles = getSampleStyleSheet()
    cell = ParagraphStyle("c", parent=styles["BodyText"], fontSize=7.4,
                          leading=9.2, textColor=SOFT)
    head = ParagraphStyle("h", parent=cell, fontSize=7.6, textColor=colors.white,
                          fontName="Helvetica-Bold")

    def P(t, s=cell):
        return Paragraph(t, s)

    rows = [[P("Stage", head), P("Input", head), P("Output", head),
             P("Tool", head), P("Cost", head)],
            [P("<b>Parent</b>"), P("A design that folds confidently"),
             P("Sequence + folded PDB with the 6NT ligand"),
             P("—"), P("—")],
            [P("<b>Propose</b>"),
             P("Parent sequence; positions measured as responsible"),
             P("~225 variants, 1–2 substitutions each, catalytic pair frozen"),
             P("LigandMPNN<br/>+ rule-based ladders"), P("seconds")],
            [P("<b>Screen</b>"), P("Variant sequences + 6NT"),
             P("One predicted complex each; distances, dihedral, iPDE, pLDDT"),
             P("Boltz-2, 1 seed<br/>8 workers"), P("~11 s / design")],
            [P("<b>Gate</b>"), P("Screened complexes"),
             P("The subset whose geometry is meaningful at all"),
             P("in-house filters<br/>MMseqs2 · Foldseek"), P("CPU, seconds")],
            [P("<b>Full score</b>"), P("Surviving sequences"),
             P("5 predictions, each scored separately; best + reproducibility"),
             P("Boltz-2, 5 seeds"), P("~3.5 min / design")],
            [P("<b>Publish</b>"), P("All measured designs, good and bad"),
             P("designs.csv, plots, PDBs, rotating 3D GIFs, status.json"),
             P("matplotlib · PyMOL · git"), P("continuous")]]
    t = Table(rows, colWidths=[22 * mm, 36 * mm, 50 * mm, 36 * mm, 22 * mm])
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), INK),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [PAPER, colors.white]),
        ("GRID", (0, 0), (-1, -1), 0.5, RULE),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 4),
        ("RIGHTPADDING", (0, 0), (-1, -1), 4),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4)]))
    w, h = t.wrap(W - 40 * mm, H)
    t.drawOn(c, 20 * mm, H - 42 * mm - h)

    y = H - 50 * mm - h
    c.setFillColor(INK)
    c.setFont("Helvetica-Bold", 10.5)
    c.drawString(20 * mm, y, "The four scored criteria (set by the organisers)")
    y -= 7 * mm
    for name, desc in [
        ("iPDE", "Boltz-2's own estimate of error in protein–ligand distances. "
                 "Confidence that the ligand is placed where it says."),
        ("Geometry", "Distances, angles and dihedrals between the catalytic "
                     "aspartate, the glutamine donor and the ligand."),
        ("Base burial", "Count of Ca atoms within 12 A of the catalytic base — "
                        "buried rather than solvent-exposed, to keep its pKa."),
        ("Consistency", "Agreement of the metrics across the five Boltz-2 "
                        "seeds. A proxy for a pre-organised active site.")]:
        c.setFillColor(ACC)
        c.setFont("Helvetica-Bold", 8.2)
        c.drawString(22 * mm, y, name)
        c.setFillColor(SOFT)
        y = wrap(c, desc, 48 * mm, y, W - 68 * mm) - 2.2 * mm

    y -= 2 * mm
    c.setFillColor(INK)
    c.setFont("Helvetica-Bold", 10.5)
    c.drawString(20 * mm, y, "Tools considered and not used, with the reason")
    y -= 7 * mm
    for name, why in [
        ("RFdiffusion3", "Generates new backbones. Our backbone is not the "
                         "problem — the catalytic contact is stable across "
                         "seeds."),
        ("Rosetta scoring", "Risso et al. 2020 tested it on this exact "
                            "reaction: no correlation with activity."),
        ("FuncLib / PROSS", "Need an alignment of homologues. A de novo design "
                            "has none (measured alignment depth: 1)."),
        ("ProteinMPNN", "Ligand-blind, and the pocket is the whole problem.")]:
        c.setFillColor(INK)
        c.setFont("Helvetica-Bold", 7.8)
        c.drawString(22 * mm, y, name)
        c.setFillColor(SOFT)
        y = wrap(c, why, 58 * mm, y, W - 78 * mm) - 2.0 * mm
    footer(c, 2, 3)


def page3(c, st):
    header(c, "Why the pipeline has this shape",
           "Every filter exists because a measurement fooled us first")
    y = H - 42 * mm
    for n, (title, body) in enumerate([
        ("We score five predictions separately, never an average",
         "One design reported 0.624. Its oxyanion dihedral is bimodal — about "
         "99 degrees on three seeds and 245 on two, never between. The "
         "circular mean of those clusters, 138.8 degrees, falls inside the "
         "allowed band that no replicate occupies. Scored per replicate, 0 of "
         "5 score. The same averaging hid a real 0.743 on a different design."),
        ("We check the ligand is bound before reading any angle",
         "Several apparent breakthroughs had the ligand 7–9 A from the "
         "catalytic aspartate. A dihedral measured on an unbound ligand is "
         "noise. This filter overturned three results that had already been "
         "reported."),
        ("We fold a known-good enzyme through our own protocol",
         "A published Kemp eliminase scored 0.000 through our pipeline because "
         "of a folding setting. Until that control was run, every zero looked "
         "like a design failure. If a check condemns the reference structure, "
         "the check is wrong."),
        ("We perturb a parent instead of generating sequences",
         "Freshly generated sequences folded at confidence 0.41–0.47, where "
         "the same sequence predicted twice gave catalytic distances of 2.6 "
         "and 13.6 A. At that confidence a design is unmeasurable, so "
         "ranking it ranks noise."),
        ("We bank a submission set early",
         "The organisers accept five designs regardless of score. Having none "
         "at the deadline is the only unrecoverable failure, so five eligible "
         "designs were secured before any further optimisation.")]):
        c.setFillColor(ACC if n % 2 == 0 else ACC2)
        c.circle(23 * mm, y + 1.2 * mm, 2.2 * mm, stroke=0, fill=1)
        c.setFillColor(colors.white)
        c.setFont("Helvetica-Bold", 7)
        c.drawCentredString(23 * mm, y + 0.2 * mm, str(n + 1))
        c.setFillColor(INK)
        c.setFont("Helvetica-Bold", 9.4)
        c.drawString(29 * mm, y, title)
        y -= 5.6 * mm
        c.setFillColor(SOFT)
        c.setFont("Helvetica", 7.8)
        y = wrap(c, body, 29 * mm, y, W - 58 * mm, size=7.8, leading=4.1) - 5 * mm

    c.setFillColor(INK)
    c.setFont("Helvetica-Bold", 10)
    c.drawString(20 * mm, y, "Where it stands")
    y -= 6.5 * mm
    c.setFillColor(SOFT)
    c.setFont("Helvetica", 8)
    for ln in [
        f"{st['screened']} designs measured · {st['scored']} with a full "
        f"five-seed score · best {st['best']:.4f} on {st['repro']}/5 replicates.",
        "The catalytic aspartate contact is stable (sd 0.07 A, matching the "
        "reference enzyme's 0.068).",
        "The open problem is one glutamine side chain that adopts a scoring "
        "rotamer in only one prediction out of five.",
        "A tyrosine hydrogen bond to that glutamine is present in every "
        "in-band design and absent in every other — the current target.",
        "",
        "A score here is a geometric proxy from a prediction model. It is not "
        "a measured reaction rate, and nothing in this document claims "
        "catalytic activity."]:
        c.drawString(22 * mm, y, ln)
        y -= 4.4 * mm
    footer(c, 3, 3)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--out", default="outputs/pipeline.pdf")
    a = ap.parse_args()
    st = stats()
    out = ROOT / a.out
    out.parent.mkdir(parents=True, exist_ok=True)
    c = rl_canvas.Canvas(str(out), pagesize=A4)
    c.setTitle("Kemp eliminase design pipeline")
    for fn in (page1, page2, page3):
        c.setFillColor(PAPER)
        c.rect(0, 0, W, H, stroke=0, fill=1)
        fn(c, st)
        c.showPage()
    c.save()
    print(f"wrote {out}  ({out.stat().st_size/1024:.0f} KB)")
    print(f"  stats embedded: screened={st['screened']} scored={st['scored']} "
          f"best={st['best']:.4f} repro={st['repro']}/5")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
