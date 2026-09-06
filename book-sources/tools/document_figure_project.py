#!/usr/bin/env python3
"""Refresh source guides and the audit index without rewriting historical logs."""
from collections import defaultdict
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
STATES = {
    'accepted': 'Accepted reconstruction, integrated into the book and independent chapters.',
    'retained-tikz': 'Previously accepted TikZ drawing retained in the reading editions.',
    'revision-requested': 'Revised candidate awaiting review; the reading editions retain the previous figure.',
    'keep-original': 'Original retained at the author’s request; the rejected proposal is archived.',
    'removed': 'Removed from the book at the author’s request.',
    'merged': 'Merged into the preceding figure; no separate comparison page.',
}


def heading(item):
    if item['chapter'] == 'cover':
        return 'Book cover'
    if not item['figure_number']:
        return 'Unnumbered illustration: ' + item['figure_slug'].replace('-', ' ')
    return 'Figure ' + item['figure_number'] + (' (continued)' if item['continued'] else '')


def previous(item):
    number = item.get('review_number')
    return str(number) if number else 'Unnumbered'


def main():
    items = json.loads((ROOT / 'docs/figures/figure-project.json').read_text())
    index = {i['id']: i for i in json.loads((ROOT / 'docs/figures/figure-index.json').read_text())}
    groups = defaultdict(list)
    for item in items:
        groups[item['code_directory']].append(item)
    for directory, figures in groups.items():
        text = '# ' + '; '.join(heading(i) for i in figures) + '\n\n'
        text += figures[0]['chapter_title'] + '.\n\n'
        for item in figures:
            state = item['review_state']
            text += '## ' + heading(item) + '\n\n' + STATES[state] + '\n\n'
            if item.get('review_number') != item['figure_number']:
                text += 'Previous audit identifier: **' + previous(item) + '**.\n\n'
            caption_label = 'Last published caption' if state in ('removed', 'merged') else 'Exact current book caption'
            text += caption_label + ' (LaTeX):\n\n```tex\n' + item['caption_latex'].strip() + '\n```\n\n'
            if state not in ('keep-original', 'removed', 'merged'):
                text += item['notes'] + '\n\n'
            if item['id'] in index:
                text += 'Comparison PDF page: **' + str(index[item['id']]['comparison_pdf_page']) + '**. '
            else:
                text += 'Omitted from the current comparison PDF. '
            text += 'Stable identifier: `' + item['id'] + '`.\n\n'
            if state not in ('keep-original', 'removed', 'merged'):
                text += 'Rebuild from the repository root:\n\n```sh\n'
                text += 'build/figure-runtime/bin/python tools/regenerate_figures.py --id ' + item['id'] + '\n```\n\n'
        text += 'Matching asset directory: `' + figures[0]['asset_directory'] + '`. '
        text += '`context.tex` records the mathematical context at reconstruction time. '
        text += '`original/` preserves historical assets and `original.pdf` assembles the previous figure. '
        text += '`proposed.pdf` is used by the reading editions only when the manifest state is `accepted`.\n'
        if figures[0]['review_state'] not in ('keep-original', 'removed', 'merged'):
            if figures[0]['engine'] == 'python':
                text += '\n`generate.py` specifies the numerical experiment; `provenance.json` records inputs, parameters, random seeds and checks.\n'
            else:
                text += '\nThe TikZ sources and `proposed.tex` specify the editable drawing and its assembly.\n'
        (ROOT / directory / 'README.md').write_text(text)

    catalog = '# Figure audit index\n\n'
    if index:
        catalog += f'Open [the separate comparison PDF](figure-comparisons.pdf). {len(index)} revised entries remain for review. '
        catalog += 'Headings use exact current book numbers and captions; the previous audit number is retained when numbering changed. '
        catalog += 'Accepted replacements, explicitly retained originals, and removed entries are omitted.\n\n'
        catalog += '| Current figure | Previous audit number | Chapter | Comparison PDF page | Source folder |\n| --- | --- | --- | ---: | --- |\n'
    else:
        catalog += 'All reviewed reconstructions are integrated into the reading editions. No comparisons remain pending. '
        catalog += 'The [previous 43-entry audit](archive/2026-09-06-reviewed-figures/figure-comparisons.pdf) and its index preserve the review history. '
        catalog += 'The current [decision ledger](author-decisions.json) records all current and previous figure numbers.\n'
    for item in items:
        if item['id'] not in index:
            continue
        old = '7.9 + 7.10' if item['id'] == 'denoising--filtering-optimal-1d' else previous(item)
        catalog += '| ' + heading(item) + ' | ' + old + ' | ' + item['chapter_title'] + ' | ' + str(index[item['id']]['comparison_pdf_page']) + ' | [' + item['figure_slug'] + '](../../' + item['code_directory'] + '/README.md) |\n'
    (ROOT / 'docs/figures/figure-index.md').write_text(catalog)
    decisions = [dict(id=i['id'], review_number=i.get('review_number'), current_number=i['figure_number'] if i['review_state'] not in ('removed','merged') else None, decision=i['review_state'], in_book=i['in_book'], comparison_pdf_page=index.get(i['id'],{}).get('comparison_pdf_page')) for i in items]
    (ROOT / 'docs/figures/author-decisions.json').write_text(json.dumps(decisions, indent=2)+'\n')
    print(f'Documented {len(items)} inventory entries in {len(groups)} source directories; {len(index)} comparisons.')


if __name__ == '__main__':
    main()
