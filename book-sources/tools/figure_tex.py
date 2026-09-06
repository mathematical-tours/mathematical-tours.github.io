"""TeX escaping and stabilized book-reference parsing."""
import re
import unicodedata

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
                    "²": r"\textsuperscript{2}", "³": r"\textsuperscript{3}", "ᵀ": r"\textsuperscript{T}", "ℓ": r"\ensuremath{\ell}", "√": r"\ensuremath{\surd}"}
    for char, command in zip("αβγδεθλμνπρστφωΣΦΩ", (
            "alpha beta gamma delta epsilon theta lambda mu nu pi rho sigma tau phi omega Sigma Phi Omega").split()):
        replacements[char] = "\\ensuremath{\\" + command + "}"
    return "".join(replacements.get(char, char) for char in unicodedata.normalize("NFC", value))



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
    chapter_aux = ([build / "book/chapters-tex" / f"{name}.aux" for name in ["abstract", *chapters]]
                   if chapters is not None else sorted((build / "book/chapters-tex").glob("*.aux")))
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
    directory = build / "figure-review"
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "book-external.aux").write_text("\n".join(lines) + "\n")
    return labels
