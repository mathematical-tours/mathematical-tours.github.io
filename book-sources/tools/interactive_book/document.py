#!/usr/bin/env python3
"""Write the publication audit guide from the generated book and checks."""
from pathlib import Path
import argparse
import collections
import datetime
import json
import re
import convert

ROOT = Path(__file__).resolve().parents[2]
WEB = ROOT / 'myst'
OUT = ROOT / 'docs/interactive-book'


def cell(text):
    return text.replace('|', '&#124;').replace('\n', ' ')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--base-url', default='http://127.0.0.1:8790',
                        help='Running book URL used by the figure audit links.')
    args = parser.parse_args()
    base = args.base_url.rstrip('/')
    figures = json.loads((WEB / 'live/figures.json').read_text())
    conversion = json.loads((WEB / 'conversion-report.json').read_text())
    metadata = json.loads((WEB / 'reference-metadata.json').read_text())
    convert.LABELS = metadata['labels']
    convert.OWNERS = metadata['owners']
    convert.BIB = metadata['bibliography']
    anchors = {}
    for label, fields in convert.LABELS.items():
        if len(fields) > 3 and fields[3]:
            key = (convert.OWNERS[label], fields[3])
            if key in anchors:
                convert.ALIASES[label] = anchors[key]
            else:
                anchors[key] = label
    renderer = convert.Renderer(figures)
    groups = collections.defaultdict(list)
    states = {}
    for figure in figures:
        groups[figure['chapter']].append(figure)
        if figure['interaction']:
            states[figure['id']] = json.loads(
                (WEB / 'live/data' / (figure['id'] + '.json')).read_text())

    lines = ['# Figure audit guide\n\n',
             'Generated from the current publication captions and experiment data. '
             f'Links open [the interactive book]({base}/); '
             'use the **Book figure** tab to compare with the approved PDF illustration. '
             'Numbers are those of the main PDF, including its continued figures.\n\n',
             'Sliders select the computed values listed below. Menus select actual '
             'inputs or methods; experiment notes explain numerical conventions. '
             'Each experiment adds an exploration prompt, suggested comparisons, and '
             'a **What to look for** explanation. **Open figure** preserves selected '
             'parameters and the comparison tab in a shareable URL. Focus a plot to '
             'inspect values with the arrow keys, or expand its numerical table. '
             'The cover features a credited painting detail by Louis Peyré; explicitly retained mathematical illustrations remain unchanged.\n\n',
             '| Chapter | Interactive | Retained | Computed states |\n',
             '| --- | ---: | ---: | ---: |\n']
    for chapter in conversion['chapters']:
        items = groups[chapter['chapter']]
        count = sum(bool(f['interaction']) for f in items)
        total = sum(len(states[f['id']]['states']) for f in items if f['interaction'])
        title = f"{chapter['number']}. {chapter['title']}"
        lines.append(f"| [{title}](#{chapter['chapter']}) | {count} | {len(items)-count} | {total} |\n")

    for chapter in conversion['chapters']:
        name = chapter['chapter']
        lines.append(f'\n<a id="{name}"></a>\n\n## {chapter["number"]}. {chapter["title"]}\n\n')
        lines.extend(['| Publication figure and exact caption | Controls and available values | States |\n',
                      '| --- | --- | ---: |\n'])
        convert.PAGE = name
        for figure in groups[name]:
            label = 'Figure ' + figure['figure_number'] if figure['figure_number'] else 'Unnumbered illustration'
            caption = renderer.inline_blocks(convert.pandoc(figure['caption_source']))
            # Captions occasionally reference another equation or figure. Resolve
            # those links into the same running edition as the audit links.
            caption = re.sub(r'\]\(([^():]*\.md)(#[^)]*)?\)',
                             lambda m: '](' + base + '/' + Path(m[1]).stem + '/' + (m[2] or '') + ')', caption)
            caption = re.sub(r'\]\((#[^)]*)\)',
                             lambda m: '](' + base + '/' + name + '/' + m[1] + ')', caption)
            target = base + '/live/panels/' + figure['id'] + '.html'
            title = f'**[{label}]({target})**<br>{cell(caption)}'
            if figure['interaction']:
                experiment = states[figure['id']]
                controls = '<br>'.join('**' + cell(c['label']) + ':** ' +
                                       ', '.join(cell(str(v)) for v in c['values'])
                                       for c in experiment['controls'])
                count = str(len(experiment['states']))
            else:
                controls = ('Retained approved TikZ display.' if figure['review_state'] == 'retained-tikz'
                            or figure['engine'] == 'tikz' else 'Original retained at the author’s request.')
                count = '—'
            lines.append(f'| {title} | {controls} | {count} |\n')
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'figures.md').write_text(''.join(lines))

    # Publish only completed, successful audits. Intermediate or failed reports
    # remain in build/ and never become a release claim.
    validation_path = ROOT / 'build/interactive-book/validation.json'
    browser_path = ROOT / 'build/interactive-book/browser/report.json'
    reading_path = ROOT / 'build/interactive-book/browser/reading-report.json'
    pdf_path = ROOT / 'build/interactive-book/pdf-links.json'
    if all(p.exists() for p in [validation_path, browser_path, reading_path, pdf_path]):
        validation = json.loads(validation_path.read_text())
        browser = json.loads(browser_path.read_text())
        reading = json.loads(reading_path.read_text())
        pdf_links = json.loads(pdf_path.read_text())
        if (validation['errors'] or browser['issues'] or browser['external_requests']
                or reading['issues'] or reading['external_requests']
                or pdf_links['errors']
                or browser['passed'] != len(figures)
                or browser['chapters_checked'] != len(conversion['chapters'])
                or not validation.get('content_sha256')
                or validation['content_sha256'] != browser.get('content_sha256')
                or validation['content_sha256'] != reading.get('content_sha256')
                or validation['content_sha256'] != pdf_links.get('content_sha256')):
            raise SystemExit('Figure guide written; release summary needs successful validation and browser checks.')
        summary = {k: v for k, v in validation.items() if k != 'interactions'}
        summary['browser'] = {k: v for k, v in browser.items() if k != 'results'}
        summary['reading'] = reading
        summary['pdf_links'] = pdf_links
        summary['validation_date'] = datetime.date.today().isoformat()
        summary['execution'] = 'Static site with explicitly computed parameter grids.'
        (OUT / 'validation.json').write_text(json.dumps(summary, indent=2, ensure_ascii=False) + '\n')
    print(f'Wrote the exact-number guide for {len(figures)} figure displays.')


if __name__ == '__main__':
    main()
