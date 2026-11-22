"""Render KEMP_FINAL_REPORT.md to PDF and DOCX.

    python scripts/render_report.py

One Markdown source, two renderers, so the three documents cannot disagree.

Deliberately a small, predictable Markdown subset - headings, tables, bullets,
numbered lists, fenced code, bold/inline-code spans, horizontal rules. No
Unicode sub/superscripts anywhere: ReportLab's built-in fonts have no glyphs
for them and they render as solid black boxes.
"""
from __future__ import annotations
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "published/report/KEMP_FINAL_REPORT.md"
OUT = ROOT / "published/report"

INLINE = re.compile(r"\*\*(.+?)\*\*|`(.+?)`|\*(.+?)\*")


def parse(md: str):
    """Markdown -> a flat list of (kind, payload) blocks."""
    blocks, lines, i = [], md.splitlines(), 0
    while i < len(lines):
        ln = lines[i]
        if ln.startswith("```"):
            buf = []
            i += 1
            while i < len(lines) and not lines[i].startswith("```"):
                buf.append(lines[i]); i += 1
            blocks.append(("code", buf)); i += 1; continue
        if ln.startswith("|") and i + 1 < len(lines) and set(
                lines[i + 1].replace("|", "").replace(" ", "")) <= set("-:"):
            hdr = [c.strip() for c in ln.strip("|").split("|")]
            i += 2
            body = []
            while i < len(lines) and lines[i].startswith("|"):
                body.append([c.strip() for c in lines[i].strip("|").split("|")])
                i += 1
            blocks.append(("table", (hdr, body))); continue
        if ln.startswith("#"):
            lv = len(ln) - len(ln.lstrip("#"))
            blocks.append((f"h{min(lv,4)}", ln.lstrip("# ").strip())); i += 1; continue
        if ln.strip() in ("---", "***", "___"):
            blocks.append(("rule", None)); i += 1; continue
        if re.match(r"^\s*[-*]\s+", ln):
            blocks.append(("li", re.sub(r"^\s*[-*]\s+", "", ln))); i += 1; continue
        if re.match(r"^\s*\d+\.\s+", ln):
            blocks.append(("oli", ln.strip())); i += 1; continue
        if ln.strip():
            buf = [ln.strip()]
            i += 1
            while (i < len(lines) and lines[i].strip()
                   and not lines[i].startswith(("#", "|", "```", "- ", "* "))
                   and not re.match(r"^\s*\d+\.\s+", lines[i])
                   and lines[i].strip() not in ("---", "***")):
                buf.append(lines[i].strip()); i += 1
            blocks.append(("p", " ".join(buf))); continue
        i += 1
    return blocks


# ------------------------------------------------------------------- PDF
def to_pdf(blocks, path: Path) -> None:
    from reportlab.lib import colors
    from reportlab.lib.enums import TA_LEFT
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer,
                                    Table, TableStyle, HRFlowable, Preformatted)

    ss = getSampleStyleSheet()
    body = ParagraphStyle("body", parent=ss["BodyText"], fontSize=9,
                          leading=12.6, spaceAfter=5, alignment=TA_LEFT)
    h1 = ParagraphStyle("h1", parent=ss["Heading1"], fontSize=18, leading=22,
                        spaceBefore=4, spaceAfter=9,
                        textColor=colors.HexColor("#15303f"))
    h2 = ParagraphStyle("h2", parent=ss["Heading2"], fontSize=13.5, leading=17,
                        spaceBefore=13, spaceAfter=6,
                        textColor=colors.HexColor("#1d566e"))
    h3 = ParagraphStyle("h3", parent=ss["Heading3"], fontSize=10.8, leading=14,
                        spaceBefore=9, spaceAfter=4,
                        textColor=colors.HexColor("#2b6f86"))
    cell = ParagraphStyle("cell", parent=body, fontSize=7.2, leading=9.0,
                          spaceAfter=0)
    cellh = ParagraphStyle("cellh", parent=cell, fontName="Helvetica-Bold",
                           textColor=colors.white)
    mono = ParagraphStyle("mono", parent=body, fontName="Courier", fontSize=8,
                          leading=10.5, leftIndent=8)

    def rich(t: str) -> str:
        """Markdown inline spans -> ReportLab markup, escaping first."""
        t = (t.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))
        t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
        t = re.sub(r"`(.+?)`", r'<font face="Courier">\1</font>', t)
        t = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"<i>\1</i>", t)
        return t

    story, first = [], True
    for kind, payload in blocks:
        if kind == "h1":
            story.append(Paragraph(rich(payload), h1)); first = False
        elif kind == "h2":
            story.append(Paragraph(rich(payload), h2))
        elif kind in ("h3", "h4"):
            story.append(Paragraph(rich(payload), h3))
        elif kind == "p":
            story.append(Paragraph(rich(payload), body))
        elif kind == "li":
            story.append(Paragraph("&bull;&nbsp;" + rich(payload), body))
        elif kind == "oli":
            story.append(Paragraph(rich(payload), body))
        elif kind == "code":
            story.append(Preformatted("\n".join(payload), mono))
            story.append(Spacer(1, 4))
        elif kind == "rule":
            if not first:
                story.append(Spacer(1, 3))
                story.append(HRFlowable(width="100%", thickness=0.6,
                                        color=colors.HexColor("#c3d3da")))
                story.append(Spacer(1, 5))
        elif kind == "table":
            hdr, rows = payload
            n = len(hdr)
            avail = A4[0] - 24 * mm
            data = [[Paragraph(rich(c), cellh) for c in hdr]]
            for r in rows:
                r = (r + [""] * n)[:n]
                data.append([Paragraph(rich(c), cell) for c in r])
            # widen the first column, which holds the design names
            if n > 4:
                w = [avail * 0.19] + [avail * 0.81 / (n - 1)] * (n - 1)
            else:
                w = [avail / n] * n
            t = Table(data, colWidths=w, repeatRows=1, hAlign="LEFT")
            t.setStyle(TableStyle([
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#2b6f86")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1),
                 [colors.white, colors.HexColor("#eef4f7")]),
                ("GRID", (0, 0), (-1, -1), 0.3, colors.HexColor("#b8ccd4")),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 3),
                ("RIGHTPADDING", (0, 0), (-1, -1), 3),
                ("TOPPADDING", (0, 0), (-1, -1), 2.5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 2.5),
            ]))
            story.append(t); story.append(Spacer(1, 7))

    doc = SimpleDocTemplate(
        str(path), pagesize=A4, leftMargin=12 * mm, rightMargin=12 * mm,
        topMargin=13 * mm, bottomMargin=13 * mm,
        title="Kemp eliminase design - final report",
        author="DBSSE Protein Design Hackathon")

    def page(canv, _doc):
        canv.saveState()
        canv.setFont("Helvetica", 7)
        canv.setFillColor(colors.HexColor("#7d969f"))
        canv.drawString(12 * mm, 7 * mm,
                        "Kemp eliminase design - final report")
        canv.drawRightString(A4[0] - 12 * mm, 7 * mm, f"page {canv.getPageNumber()}")
        canv.restoreState()

    doc.build(story, onFirstPage=page, onLaterPages=page)


# ------------------------------------------------------------------- DOCX
def to_docx(blocks, path: Path) -> None:
    from docx import Document
    from docx.shared import Pt, RGBColor, Inches
    from docx.enum.text import WD_ALIGN_PARAGRAPH

    doc = Document()
    for s in doc.sections:
        s.left_margin = s.right_margin = Inches(0.6)
        s.top_margin = s.bottom_margin = Inches(0.6)
    normal = doc.styles["Normal"]
    normal.font.name = "Calibri"
    normal.font.size = Pt(9.5)

    def runs(par, text, size=None, base_bold=False):
        """Emit runs honouring **bold**, `code` and *italic*."""
        pos = 0
        for m in INLINE.finditer(text):
            if m.start() > pos:
                r = par.add_run(text[pos:m.start()])
                r.bold = base_bold
                if size:
                    r.font.size = size
            b, c, i = m.group(1), m.group(2), m.group(3)
            r = par.add_run(b or c or i)
            r.bold = base_bold or bool(b)
            r.italic = bool(i)
            if c:
                r.font.name = "Consolas"
            if size:
                r.font.size = size
            pos = m.end()
        if pos < len(text):
            r = par.add_run(text[pos:])
            r.bold = base_bold
            if size:
                r.font.size = size

    for kind, payload in blocks:
        if kind == "h1":
            p = doc.add_heading(level=0)
            runs(p, payload)
        elif kind == "h2":
            p = doc.add_heading(level=1)
            runs(p, payload)
        elif kind in ("h3", "h4"):
            p = doc.add_heading(level=2)
            runs(p, payload)
        elif kind == "p":
            runs(doc.add_paragraph(), payload)
        elif kind == "li":
            runs(doc.add_paragraph(style="List Bullet"), payload)
        elif kind == "oli":
            runs(doc.add_paragraph(style="List Number"),
                 re.sub(r"^\s*\d+\.\s+", "", payload))
        elif kind == "code":
            p = doc.add_paragraph()
            r = p.add_run("\n".join(payload))
            r.font.name = "Consolas"
            r.font.size = Pt(8.5)
        elif kind == "rule":
            doc.add_paragraph()
        elif kind == "table":
            hdr, rows = payload
            n = len(hdr)
            t = doc.add_table(rows=1, cols=n)
            t.style = "Light Grid Accent 1"
            t.autofit = True
            for j, c in enumerate(hdr):
                cellp = t.rows[0].cells[j].paragraphs[0]
                runs(cellp, c, size=Pt(7.5), base_bold=True)
            for r in rows:
                r = (r + [""] * n)[:n]
                cells = t.add_row().cells
                for j, c in enumerate(r):
                    runs(cells[j].paragraphs[0], c, size=Pt(7.5))
            doc.add_paragraph()
    doc.save(str(path))


if __name__ == "__main__":
    blocks = parse(SRC.read_text())
    to_pdf(blocks, OUT / "KEMP_FINAL_REPORT.pdf")
    to_docx(blocks, OUT / "KEMP_FINAL_REPORT.docx")
    import subprocess
    for f in ("KEMP_FINAL_REPORT.pdf", "KEMP_FINAL_REPORT.docx"):
        p = OUT / f
        print(f"  {f:28} {p.stat().st_size/1024:7.0f} KB")
    print(f"blocks: {len(blocks)}  tables: {sum(1 for k,_ in blocks if k=='table')}")
