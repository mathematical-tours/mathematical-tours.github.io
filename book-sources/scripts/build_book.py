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
from tikz_reviews import build_reviews, prepare_comparisons

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "build"
MAIN = "FundationsDataScience"
DIAGNOSTIC = re.compile(r"^(?:LaTeX(?: Font)? Warning:|Package .+ Warning:|pdfTeX warning|Overfull \\[hv]box|Underfull \\[hv]box|Missing character:|Warning--)")

def active_text(path):
    text = re.sub(r"(?<!\\)%[^\n]*", "", path.read_text(encoding="utf-8"))
    return re.sub(r"\\iffalse\b.*?\\fi\b", "", text, flags=re.S)

def chapter_list():
    return [n for n in re.findall(r"\\include\{chapters/([^}]+)\}", active_text(ROOT / f"{MAIN}.tex")) if n != "abstract"]

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
    (directory / "chapters").mkdir(exist_ok=True)
    command = ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "-file-line-error", "-recorder",
               f"-jobname={name}", f"-output-directory={directory.relative_to(ROOT)}", str(source)]
    run(command, directory / "latex-1.txt")
    aux_paths = [directory / f"{name}.aux", *(directory / "chapters").glob("*.aux")]
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
        path = BUILD / "book" / "chapters" / f"{chapter}.aux"
        if chapter != name and path.exists():
            labels.extend(s for s in path.read_text().splitlines() if s.startswith("\\newlabel{"))
    external = directory / "book-external.aux"
    external.write_text("\n".join(labels) + "\n")
    driver = directory / "driver.tex"
    definitions = [rf"\def\StandaloneChapter{{{name}}}", rf"\def\StandaloneNumber{{{number}}}",
                   rf"\def\BookExternalReferences{{{external.with_suffix('').relative_to(ROOT)}}}"]
    if not has_citations(ROOT / "chapters" / f"{name}.tex"):
        definitions.append(r"\def\SkipBibliography{1}")
    driver.write_text("\n".join(definitions) + f"\n\\input{{{MAIN}.tex}}\n")
    return compile_document(name, driver.relative_to(ROOT), directory)

def build_figures():
    directory = BUILD / "figures"
    directory.mkdir(parents=True, exist_ok=True)
    run([sys.executable, "scripts/variance_stabilization.py"], directory / "variance-data.txt")
    for relative in ("figures/ml/classif/losses-corrected", "figures/denoising/variance-stabilization-corrected",
                     "figures/wavelets/wavelets-2d-support-corrected"):
        source = ROOT / (relative + ".tex")
        destination = source.with_suffix(".pdf")
        dependencies = [source]
        if source.with_suffix(".dat").exists():
            dependencies.append(source.with_suffix(".dat"))
        if destination.exists() and destination.stat().st_mtime >= max(p.stat().st_mtime for p in dependencies):
            continue
        run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
             f"-output-directory={directory.relative_to(ROOT)}", str(source)], directory / f"{source.stem}.txt")
        log = (directory / f"{source.stem}.log").read_text(errors="replace")
        if any(DIAGNOSTIC.match(line) for line in log.splitlines()):
            raise RuntimeError(f"Figure diagnostics in {directory / (source.stem + '.log')}")
        shutil.copy2(directory / destination.name, destination)

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
    build_figures()
    build_reviews(ROOT, BUILD, chapter_list(), run, DIAGNOSTIC, args.jobs)
    print("Building the complete book...", flush=True)
    book = compile_document(MAIN, ROOT / f"{MAIN}.tex", BUILD / "book")
    print(f"Book: {book['pages']} pages, {len(book['diagnostics'])} diagnostics", flush=True)
    prepare_comparisons(ROOT, BUILD, chapter_list())
    print("Building the separate figure comparison volume...", flush=True)
    comparisons = compile_document("figure-comparisons", ROOT / "figure-processing/figure-comparisons.tex",
                                   BUILD / "figure-processing")
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
    (ROOT / "figure-processing").mkdir(exist_ok=True)
    for report in reports:
        destination = ROOT / f"{MAIN}.pdf" if report["name"] == MAIN else ROOT / "chapters-pdf" / f"{report['name']}.pdf"
        if report["name"] == "figure-comparisons":
            destination = ROOT / "figure-processing/figure-comparisons.pdf"
        shutil.copy2(report["pdf"], destination)
    shutil.copy2(BUILD / "figure-processing/figure-index.json", ROOT / "figure-processing/figure-index.json")
    print(f"Published {len(reports)} PDFs. Report: {report_path}")
    return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, UnicodeError) as error:
        print(error, file=sys.stderr)
        raise SystemExit(1)
