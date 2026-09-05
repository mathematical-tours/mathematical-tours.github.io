"""Build editable figures and a comparison volume keyed to the compiled book."""
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import re
import shutil
import unicodedata


def entries(root):
    result = []
    for manifest in sorted((root / "figures/tikz").glob("*.json")):
        result.extend(json.loads(manifest.read_text()))
    ids = [entry["id"] for entry in result]
    if len(ids) != len(set(ids)):
        raise RuntimeError("Duplicate figure review identifiers")
    for entry in result:
        for field in ("id", "chapter", "original", "tex", "title", "context", "notes", "confidence", "book_label", "book_panel"):
            if field not in entry:
                raise RuntimeError(f"Missing {field} in figure review {entry}")
        for path in (root / "figures" / entry["original"], root / entry["tex"]):
            if not path.is_file():
                raise RuntimeError(f"Missing figure review source: {path}")
    return result


def tex_text(value):
    replacements = {"&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#",
                    "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
                    "^": r"\textasciicircum{}", "\\": r"\textbackslash{}",
                    "–": "--", "—": "---", "’": "'", "‘": "`", "“": "``", "”": "''",
                    "≤": r"\ensuremath{\leq}", "≥": r"\ensuremath{\geq}",
                    "±": r"\ensuremath{\pm}", "×": r"\ensuremath{\times}",
                    "→": r"\ensuremath{\to}", "∞": r"\ensuremath{\infty}",
                    "∈": r"\ensuremath{\in}", "∥": r"\ensuremath{\Vert}",
                    "〈": r"\ensuremath{\langle}", "〉": r"\ensuremath{\rangle}", "−": "-",
                    "₁": r"\textsubscript{1}", "₂": r"\textsubscript{2}",
                    "²": r"\textsuperscript{2}", "³": r"\textsuperscript{3}"}
    for char, command in zip("αβγδεθλμνπρστφωΣΦΩ", (
            "alpha beta gamma delta epsilon theta lambda mu nu pi rho sigma tau phi omega Sigma Phi Omega").split()):
        replacements[char] = "\\ensuremath{\\" + command + "}"
    return "".join(replacements.get(char, char) for char in unicodedata.normalize("NFC", value))


def context_text(value):
    """Render book-label tokens as numbered, clickable mathematical references."""
    pattern = r"\b(?:eq|fig|prop|thm|sec|lem)-[A-Za-z0-9_:-]+"
    parts = re.split("(" + pattern + ")", value)
    return "".join((r"\eqref{" if part.startswith("eq-") else r"\ref{") + part + "}"
                   if re.fullmatch(pattern, part) else tex_text(part) for part in parts)


def build_reviews(root, build, chapters, run, diagnostic, jobs=4):
    figures = entries(root)
    directory = build / "tikz-figures"
    directory.mkdir(parents=True, exist_ok=True)
    preamble = root / "figures/tikz/tikz-preamble.tex"

    def compile_figure(entry):
        source = root / entry["tex"]
        destination = source.with_suffix(".pdf")
        if destination.exists() and destination.stat().st_mtime >= max(source.stat().st_mtime, preamble.stat().st_mtime):
            return
        run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "-file-line-error",
             f"-output-directory={directory.relative_to(root)}", str(source)], directory / f"{source.stem}.txt")
        log = (directory / f"{source.stem}.log").read_text(errors="replace")
        problems = [line for line in log.splitlines() if diagnostic.match(line)]
        if problems:
            raise RuntimeError(f"Figure diagnostics in {source.name}:\n" + "\n".join(problems))
        shutil.copy2(directory / destination.name, destination)

    sources = {entry["tex"]: entry for entry in figures}
    with ThreadPoolExecutor(max_workers=max(1, jobs)) as pool:
        list(pool.map(compile_figure, sources.values()))

    unknown = {entry["chapter"] for entry in figures} - set(chapters)
    if unknown:
        raise RuntimeError(f"Unknown figure-review chapters: {unknown}")
    print(f"Prepared {len(figures)} TikZ reconstructions.", flush=True)


def braced(text, position):
    """Read one balanced TeX group, preserving nested math and escaped braces."""
    while position < len(text) and text[position].isspace():
        position += 1
    if position >= len(text) or text[position] != "{":
        raise RuntimeError(f"Expected a TeX group at {text[position:position + 80]!r}")
    start = position + 1
    depth = 1
    position += 1
    while position < len(text) and depth:
        if text[position] == "\\":
            position += 2
            continue
        if text[position] == "{":
            depth += 1
        elif text[position] == "}":
            depth -= 1
        position += 1
    if depth:
        raise RuntimeError("Unbalanced TeX group in book references")
    return text[start:position - 1], position


def book_labels(root, build, chapters=None):
    """Read the stabilized book's labels, including full caption math."""
    labels = {}
    lines = []
    chapter_aux = ([build / "book/chapters" / f"{name}.aux" for name in ["abstract", *chapters]]
                   if chapters is not None else sorted((build / "book/chapters").glob("*.aux")))
    for path in [build / "book/FundationsDataScience.aux", *chapter_aux]:
        for line in path.read_text().splitlines():
            if not line.startswith(r"\newlabel{"):
                continue
            label, position = braced(line, len(r"\newlabel"))
            body, _ = braced(line, position)
            fields = []
            position = 0
            while position < len(body):
                field, position = braced(body, position)
                fields.append(field)
            if label in labels and labels[label] != fields:
                raise RuntimeError(f"Conflicting book label: {label}")
            labels[label] = fields
            lines.append(line)
    directory = build / "figure-processing"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "book-external.aux").write_text("\n".join(lines) + "\n")
    return labels


def prepare_comparisons(root, build, chapters):
    """Generate review pages only after the main book's numbering has settled."""
    figures = entries(root)
    labels = book_labels(root, build, chapters)
    review_dir = root / "figures/tikz/reviews"
    review_dir.mkdir(exist_ok=True)
    index = []
    for chapter in chapters:
        selected = [entry for entry in figures if entry["chapter"] == chapter]
        for entry in selected:
            label = entry["book_label"]
            if label not in labels or not labels[label][3].startswith("figure."):
                raise RuntimeError(f"Missing numbered book figure for {entry['id']}: {label}")
        selected.sort(key=lambda e: (tuple(int(n) for n in labels[e["book_label"]][0].split(".")),
                                     int(e["book_panel"].split()[1]) if e["book_panel"].startswith("Panel ") else 0))
        lines = ["% Generated by scripts/tikz_reviews.py using the compiled book's labels."]
        if selected:
            chapter_number = labels[selected[0]["book_label"]][0].split(".")[0]
            lines.append(rf"\pdfbookmark[0]{{Chapter {chapter_number}}}{{comparison-chapter-{chapter_number}}}")
        for entry in selected:
            label = entry["book_label"]
            number, page, caption = labels[label][:3]
            args = [label, tex_text(entry["title"]), "figures/" + entry["original"],
                    str(Path(entry["tex"]).with_suffix(".pdf")), context_text(entry["context"]),
                    tex_text(("Interpretation to check: " if entry["confidence"] == "qualified" else "") + entry["notes"]),
                    entry["id"], tex_text(entry["book_panel"]), number]
            lines.append("\\ReviewFigure" + "\n".join("{" + arg + "}" for arg in args))
            index.append(dict(id=entry["id"], chapter=chapter, figure_number=number,
                              comparison_pdf_page=len(index) + 2,
                              book_page=page, book_label=label, caption_latex=caption,
                              panel=entry["book_panel"], in_book=entry.get("in_book", False),
                              reconstruction_title=entry["title"],
                              original=entry["original"], reconstruction=entry["tex"]))
        target = review_dir / f"{chapter}.tex"
        if selected:
            target.write_text("\n".join(lines) + "\n")
        elif target.exists():
            target.unlink()
    includes = [rf"\input{{figures/tikz/reviews/{chapter}}}" for chapter in chapters
                if any(e["chapter"] == chapter for e in figures)]
    (build / "figure-processing/comparison-pages.tex").write_text("\n".join(includes) + "\n")
    (build / "figure-processing/figure-index.json").write_text(json.dumps(index, indent=2, ensure_ascii=False) + "\n")
    return index
