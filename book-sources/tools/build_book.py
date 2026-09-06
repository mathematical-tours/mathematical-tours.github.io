#!/usr/bin/env python3
"""Build the book and standalone chapters from the active include list."""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from figure_tex import book_labels

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "build"
MAIN = "FundationsDataScience"
DIAGNOSTIC = re.compile(r"^(?:LaTeX(?: Font)? Warning:|Package .+ Warning:|pdfTeX warning|Overfull \\[hv]box|Underfull \\[hv]box|Missing character:|Warning--)")

def active_text(path):
    text = re.sub(r"(?<!\\)%[^\n]*", "", path.read_text(encoding="utf-8"))
    return re.sub(r"\\iffalse\b.*?\\fi\b", "", text, flags=re.S)

def chapter_list():
    return [n for n in re.findall(r"\\include\{chapters-tex/([^}]+)\}", active_text(ROOT / f"{MAIN}.tex")) if n != "abstract"]

def has_citations(path):
    text = active_text(path)
    if re.search(r"\\(?:cite\w*|nocite)\b", text):
        return True
    return any(has_citations(ROOT / (n if n.endswith(".tex") else n + ".tex"))
               for n in re.findall(r"\\input\{([^}]+)\}", text))

def run(command, console, cwd=ROOT):
    env = os.environ.copy()
    env["max_print_line"] = "1000"
    env["BIBINPUTS"] = str(ROOT) + os.pathsep + env.get("BIBINPUTS", "")
    with console.open("w") as out:
        result = subprocess.run(command, cwd=cwd, env=env, stdout=out, stderr=subprocess.STDOUT, check=False)
    if result.returncode:
        tail = "\n".join(console.read_text(errors="replace").splitlines()[-25:])
        raise RuntimeError(f"Build failed; see {console}\n{tail}")

def fingerprint(directory):
    digest = hashlib.sha256()
    for path in sorted(directory.rglob("*")):
        if path.suffix in {".aux", ".toc", ".out"}:
            digest.update(path.read_bytes())
    return digest.hexdigest()

def compile_document(name, source, directory):
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "chapters-tex").mkdir(exist_ok=True)
    command = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "-file-line-error", "-recorder",
               f"-jobname={name}", f"-output-directory={directory.relative_to(ROOT)}", str(source)]
    run(command, directory / "latex-1.txt")
    aux_paths = [directory / f"{name}.aux", *(directory / "chapters-tex").glob("*.aux")]
    if any("\\citation{" in p.read_text(errors="replace") for p in aux_paths):
        run(["bibtex", name], directory / "bibtex.txt", cwd=directory)
    previous = fingerprint(directory)
    for iteration in range(2, 7):
        run(command, directory / f"latex-{iteration}.txt")
        current = fingerprint(directory)
        if current == previous and iteration >= 3:
            break
        previous = current
    else:
        raise RuntimeError(f"References did not stabilize for {name}")
    log = (directory / f"{name}.log").read_text(errors="replace")
    blg = directory / f"{name}.blg"
    diagnostics = [s for s in (log + (blg.read_text(errors="replace") if blg.exists() else "")).splitlines() if DIAGNOSTIC.match(s)]
    pages = re.search(r"Output written on .*?\((\d+) pages?", log, re.S)
    return {"name": name, "pages": int(pages.group(1)) if pages else None, "diagnostics": diagnostics,
            "passes": iteration, "pdf": str(directory / f"{name}.pdf")}

def build_chapter(name, number):
    directory = BUILD / "chapters" / name
    directory.mkdir(parents=True, exist_ok=True)
    # Import other chapters' labels only; local labels remain authoritative.
    labels = []
    for chapter in ["abstract", *chapter_list()]:
        path = BUILD / "book" / "chapters-tex" / f"{chapter}.aux"
        if chapter != name and path.exists():
            labels.extend(s for s in path.read_text().splitlines() if s.startswith("\\newlabel{"))
    external = directory / "book-external.aux"
    # xr-hyper must recognize the already complete five-field label records;
    # otherwise \hyperref links acquire spurious compatibility fields.
    external.write_text("\\HyperFirstAtBeginDocument{}\n" + "\n".join(labels) + "\n")
    driver = directory / "driver.tex"
    definitions = [rf"\def\StandaloneChapter{{{name}}}", rf"\def\StandaloneNumber{{{number}}}",
                   rf"\def\BookExternalReferences{{{external.with_suffix('').relative_to(ROOT)}}}"]
    if not has_citations(ROOT / "chapters-tex" / f"{name}.tex"):
        definitions.append(r"\def\SkipBibliography{1}")
    driver.write_text("\n".join(definitions) + f"\n\\input{{{MAIN}.tex}}\n")
    return compile_document(name, driver.relative_to(ROOT), directory)

def build_figures(jobs=4):
    directory=BUILD/'figure-regeneration';directory.mkdir(parents=True,exist_ok=True)
    python=os.environ.get('FIGURE_PYTHON')
    if not python:
        local=BUILD/'figure-runtime/bin/python'
        python=str(local) if local.exists() else sys.executable
    run([python,'tools/regenerate_figures.py','--jobs',str(jobs)],directory/'regeneration-console.txt')

def ensure_book_references():
    """Seed figure label references from checked-in assets on a fresh checkout."""
    if not (BUILD / "figure-review/book-external.aux").exists():
        compile_document(MAIN, ROOT / f"{MAIN}.tex", BUILD / "book")
        book_labels(ROOT, BUILD, chapter_list())

def check_reading_assets():
    items=json.loads((ROOT/'docs/figures/figure-project.json').read_text())
    for item in items:
        if item.get('in_book') and item.get('review_state')!='accepted':raise RuntimeError('Only author-accepted reconstructions may be integrated: '+item['id'])
        for asset in item['published_assets']:
            path=ROOT/'figures'/asset['published']
            if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest()!=asset['sha256']:
                raise RuntimeError(f'Preserved reading asset changed or missing: {path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--book-only", action="store_true")
    parser.add_argument("--chapter", choices=chapter_list())
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--allow-warnings", action="store_true", help="Publish development PDFs despite diagnostics")
    args = parser.parse_args()
    for program in ("pdflatex", "bibtex"):
        if not shutil.which(program):
            parser.error(f"{program} is required; install TeX Live")
    check_reading_assets()
    build_figures(args.jobs)
    print("Building the complete book...", flush=True)
    book = compile_document(MAIN, ROOT / f"{MAIN}.tex", BUILD / "book")
    print(f"Book: {book['pages']} pages, {len(book['diagnostics'])} diagnostics", flush=True)
    book_labels(ROOT, BUILD, chapter_list())
    print("Building the separate figure comparison volume...", flush=True)
    from build_figure_comparator import build_comparisons
    comparisons=build_comparisons(preview=True)
    comparisons["name"]="figure-comparisons"
    reports = [book, comparisons]
    print(f"Comparisons: {comparisons['pages']} pages, {len(comparisons['diagnostics'])} diagnostics", flush=True)
    if not args.book_only:
        selected = [(n, i) for i, n in enumerate(chapter_list(), 1) if args.chapter is None or args.chapter == n]
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            futures = [pool.submit(build_chapter, name, number) for name, number in selected]
            for future in as_completed(futures):
                result = future.result()
                reports.append(result)
                print(f"{result['name']}: {result['pages']} pages, {len(result['diagnostics'])} diagnostics", flush=True)
    report_path = BUILD / "build-report.json"
    report_path.write_text(json.dumps(reports, indent=2) + "\n")
    if any(r["diagnostics"] for r in reports) and not args.allow_warnings:
        print(f"Unresolved diagnostics; PDFs have not been published. See {report_path}", file=sys.stderr)
        return 1
    (ROOT / "chapters-pdf").mkdir(exist_ok=True)
    (ROOT / "docs/figures").mkdir(parents=True, exist_ok=True)
    for report in reports:
        destination = ROOT / f"{MAIN}.pdf" if report["name"] == MAIN else ROOT / "chapters-pdf" / f"{report['name']}.pdf"
        if report["name"] == "figure-comparisons":
            destination = ROOT / "docs/figures/figure-comparisons.pdf"
        shutil.copy2(report["pdf"], destination)
    shutil.copy2(BUILD / "figure-regeneration/complete-figure-index.json", ROOT / "docs/figures/figure-index.json")
    print(f"Published {len(reports)} PDFs. Report: {report_path}")
    return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, UnicodeError) as error:
        print(error, file=sys.stderr)
        raise SystemExit(1)
