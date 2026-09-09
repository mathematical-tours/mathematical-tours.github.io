#!/usr/bin/env python3
"""Check the compiled book, registered figures/programs, and audit numbering."""
import json
from pathlib import Path
import re
from pypdf import PdfReader
from build_comparator import labels_from_aux, caption_from_source

ROOT = Path(__file__).resolve().parents[1]


def verify_exercises(book, labels):
    """Check chapter numbering and the placement of every printed solution."""
    main = (ROOT/'DAFT-EN.tex').read_text()
    chapters = re.findall(r'\\include\{(chapters-tex/[1-8]-[^}]+)\}', main)
    assert len(chapters) == 8
    destinations = [item for item in book.outline if not isinstance(item, list)]
    titles = [item.title for item in destinations]
    assert 'Solutions to the exercises' not in titles, 'Solutions still form a separate appendix.'
    page_text = [page.extract_text() for page in book.pages]
    solution_pattern = r'Solution to Exercise\s*(\d+\.\d+)\b'
    for number, text in enumerate(page_text, 1):
        assert not re.search(solution_pattern+r'\s*\d+\s*\Z', text), \
            f'A solution heading is stranded at the bottom of PDF page {number}.'
    expected_all = []
    report = []
    for number, relative in enumerate(chapters, 1):
        source = (ROOT/(relative+'.tex')).read_text()
        title = re.search(r'\\chapter\{([^}]+)\}', source)[1]
        chapter_label = re.search(r'\\label\{([^}]+)\}', source)[1]
        blocks = re.findall(r'\\begin\{exo\}.*?\\end\{exo\}', source, re.S)
        exercise_ids = [re.search(r'\\label\{([^}]+)\}', block)[1] for block in blocks]
        assert len(exercise_ids) == len(set(exercise_ids))
        expected = [f'{number}.{i}' for i in range(1, len(blocks)+1)]
        assert [labels[key][0] for key in exercise_ids] == expected
        last_input = re.search(r'\\input\{(chapters-tex/solutions/[^}]+)\}\s*\Z', source)
        assert last_input, f'{relative}: solutions must close the chapter.'
        solutions = (ROOT/(last_input[1]+'.tex')).read_text()
        assert '\\section{Solutions to the exercises}' in solutions
        assert re.findall(r'\\begin\{correction\}\{([^}]+)\}', solutions) == exercise_ids
        section_label = re.search(r'\\label\{([^}]+)\}', solutions)[1]
        assert labels[section_label][0].startswith(f'{number}.')
        position = titles.index(title)
        chapter_page = book.get_destination_page_number(destinations[position])
        next_page = book.get_destination_page_number(destinations[position+1])
        offset = chapter_page-int(labels[chapter_label][1])
        solution_page = int(labels[section_label][1])+offset
        assert chapter_page < solution_page < next_page
        assert all(int(labels[key][1])+offset < solution_page for key in exercise_ids)
        assert not re.findall(solution_pattern, '\n'.join(page_text[chapter_page:solution_page]))
        printed = re.findall(solution_pattern, '\n'.join(page_text[solution_page:next_page]))
        assert printed == expected, (relative, printed, expected)
        report.append({'chapter': number, 'exercises': len(expected),
                       'number_range': f'{expected[0]}–{expected[-1]}',
                       'chapter_pdf_page': chapter_page+1,
                       'solutions_section': labels[section_label][0],
                       'solutions_pdf_page': solution_page+1})
        expected_all.extend(expected)
    assert len(expected_all) == 88
    assert re.findall(solution_pattern, '\n'.join(page_text)) == expected_all
    assert labels['appendix-python-numerical'][0] == 'A'
    assert labels['appendix-software-resources'][0] == 'B'
    return report


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
    exercises = verify_exercises(book, labels)
    report = {'book_pages': len(book.pages), 'book_bytes': (ROOT/'DAFT-EN.pdf').stat().st_size,
              'figures': len(figures), 'programs': len(programs),
              'comparison_pages': len(audit.pages), 'tex_diagnostics': diagnostics,
              'figure_numbers_and_captions_match': True,
              'all_regenerated_figures_in_book': True,
              'printed_program_numbers': '1 through 88, each once',
              'exercises': 88, 'solutions': 88,
              'chapter_end_solutions': exercises,
              'python_appendices': ['A', 'B']}
    (ROOT/'docs/document-verification.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
