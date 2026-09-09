#!/usr/bin/env python3
"""Check the compiled book, registered figures/programs, and audit numbering."""
import json
from pathlib import Path
import re
from pypdf import PdfReader
from build_comparator import labels_from_aux, caption_from_source

ROOT = Path(__file__).resolve().parents[1]


def main():
    log = (ROOT/'DAFT-EN.log').read_text()
    diagnostics = [line for line in log.splitlines() if re.search(
        r'Warning|Overfull|Underfull|Missing character|^!', line)]
    assert not diagnostics, '\n'.join(diagnostics)
    labels = labels_from_aux(ROOT/'DAFT-EN.aux')
    figures = json.loads((ROOT/'docs/figure-manifest.json').read_text())
    programs = json.loads((ROOT/'docs/listing-manifest.json').read_text())
    comparator = json.loads((ROOT/'comparator-figures/index.json').read_text())
    actual_inputs = (ROOT/'DAFT-EN.fls').read_text()
    assert len(figures) == len(comparator) == 68
    assert len(programs) == 88
    for figure, pair in zip(figures, comparator):
        assert figure['status'] == 'regenerated'
        assert figure['output'] in actual_inputs, figure['id']
        assert (ROOT/figure['code']).is_file()
        assert (pair['number'], pair['book_page']) == labels[figure['id']]
        caption = re.sub(r'\\(eqref|ref)\{([^}]+)\}',
                         lambda m: '('+labels[m[2]][0]+')' if m[1]=='eqref'
                         else labels[m[2]][0], caption_from_source(figure))
        assert pair['caption'] == caption, figure['id']
    numbers = []
    for program in programs:
        assert program['status'] == 'python'
        assert 'listings/tex/'+program['id']+'.tex' in actual_inputs
        numbers.append(int(labels[program['id']][0]))
    assert sorted(numbers) == list(range(1, 89))
    book = PdfReader(ROOT/'DAFT-EN.pdf')
    audit = PdfReader(ROOT/'comparator-figures/comparison.pdf')
    assert len(audit.pages) == 69
    text = '\n'.join(page.extract_text() for page in book.pages)
    assert '??' not in text, 'Unresolved reference visible in the PDF.'
    assert '9782729818678' in text
    assert 'adapted English translation' in text
    assert all(len(page.extract_text().strip()) > 10 for page in book.pages)
    report = {'book_pages': len(book.pages), 'book_bytes': (ROOT/'DAFT-EN.pdf').stat().st_size,
              'figures': len(figures), 'programs': len(programs),
              'comparison_pages': len(audit.pages), 'tex_diagnostics': diagnostics,
              'figure_numbers_and_captions_match': True,
              'all_regenerated_figures_in_book': True,
              'printed_program_numbers': '1 through 88, each once'}
    (ROOT/'docs/document-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
