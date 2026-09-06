#!/usr/bin/env python3
"""Regenerate the Gaussian-maximum plot with its x-axis label below the plot.

Optional regeneration dependency: pypdf (``python -m pip install pypdf``).
Run from any directory: ``python scripts/correct_max_gaussian_std.py``.
The normal LaTeX build uses the checked-in PDF and does not run this script.

The original PDF contains an empirical plot in one embedded image, followed
by a vector-outline label. The image bytes, placement, and original drawing
commands are preserved. Only a translation around the label and a larger
page box are added. No data values are inferred or regenerated, and no fonts
are introduced: the original label already consists entirely of outlines.
"""

from hashlib import sha256
from pathlib import Path

from pypdf import PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject, NameObject, RectangleObject


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "figures/denoising/max-gaussian-std.pdf"
DESTINATION = ROOT / "figures/denoising/max-gaussian-std-corrected.pdf"
SOURCE_SHA256 = "97cb58181e7154212820624b0fb2495356ea3f46f55009d1019c512f7886c751"
LABEL_START = b"q 601 49 113.301 32 re W\nn /Gs3 gs"
LABEL_TRANSLATION = b"q\n1 0 0 1 -263.65 -90 cm\n"


def regenerate() -> None:
    source_bytes = SOURCE.read_bytes()
    if sha256(source_bytes).hexdigest() != SOURCE_SHA256:
        raise ValueError("The original figure changed; inspect its content before editing.")

    reader = PdfReader(SOURCE)
    if len(reader.pages) != 1:
        raise ValueError("Expected a one-page figure.")
    original = reader.pages[0]
    content = original.get_contents().get_data()
    if content.count(LABEL_START) != 1:
        raise ValueError("Could not identify the vector label uniquely.")
    plot_commands, label_tail = content.split(LABEL_START)
    label_commands = LABEL_START + label_tail
    if not plot_commands.rstrip().endswith(b"/Im1 Do Q Q"):
        raise ValueError("Unexpected plot drawing commands.")
    if b" Do" in label_commands or b"BT" in label_commands:
        raise ValueError("Expected only vector outlines in the label.")

    # Original label bounds: x=601..714.301, y=49..81 points.
    # Center its outlines under the plotting area, with their top at y=-9.
    # The enlarged lower page edge at -49 leaves eight points of padding.
    corrected_content = (
        plot_commands + b"\n" + LABEL_TRANSLATION + label_commands + b"\nQ\n"
    )
    writer = PdfWriter()
    page = writer.add_page(original)
    stream = DecodedStreamObject()
    stream.set_data(corrected_content)
    page[NameObject("/Contents")] = writer._add_object(stream)
    page.mediabox = RectangleObject([0, -49, 732.4748, 213])
    page.cropbox = RectangleObject([0, -49, 732.4748, 213])
    writer.write(DESTINATION)

    # Check that writing the document preserved the original plot exactly.
    saved = PdfReader(DESTINATION).pages[0]
    assert saved.get_contents().get_data() == corrected_content
    old_image = original["/Resources"]["/XObject"]["/Im1"].get_object()
    new_image = saved["/Resources"]["/XObject"]["/Im1"].get_object()
    assert old_image._data == new_image._data
    assert old_image.get_data() == new_image.get_data()
    for key in old_image:
        if key != "/ColorSpace":
            assert old_image[key] == new_image[key]
    old_space, new_space = old_image["/ColorSpace"], new_image["/ColorSpace"]
    assert old_space[0] == new_space[0] == "/ICCBased"
    assert old_space[1].get_object().get_data() == new_space[1].get_object().get_data()
    assert "/Font" not in original["/Resources"]
    assert "/Font" not in saved["/Resources"]
    assert SOURCE.read_bytes() == source_bytes
    print(f"Created {DESTINATION.relative_to(ROOT)}")
    print("Verified unchanged plot image and drawing commands; no font resources.")


if __name__ == "__main__":
    regenerate()
