#!/usr/bin/env python3
"""Check the public landing pages and local download/asset links after Jekyll builds."""
import argparse
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlsplit, unquote
import json
import sys

PAGES = ['index.html', 'book/index.html', 'daft/index.html',
         'book-basics/index.html', 'algorithms/index.html', 'slides/index.html',
         'maths-ia-course/index.html', '404.html', 'tweets/index.html',
         'tweets_main/index.html', 'tweets_oldies/index.html']

class Document(HTMLParser):
    def __init__(self, text):
        super().__init__(convert_charrefs=True)
        self.links = []
        self.ids = set()
        self.duplicate_ids = []
        self.headings = 0
        self.images_without_alt = 0
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'h1': self.headings += 1
        if tag == 'img' and 'alt' not in attrs: self.images_without_alt += 1
        if 'id' in attrs:
            if attrs['id'] in self.ids: self.duplicate_ids.append(attrs['id'])
            self.ids.add(attrs['id'])
        if tag in ('a', 'link') and attrs.get('href'): self.links.append(attrs['href'])
        if tag in ('img', 'script') and attrs.get('src'): self.links.append(attrs['src'])

def check(root):
    errors, checked = [], set()
    documents = {}
    def document(path):
        if path not in documents: documents[path] = Document(path.read_text())
        return documents[path]
    for name in PAGES:
        file = root/name
        if not file.is_file():
            errors.append(f'Missing page: {name}')
            continue
        doc = document(file)
        if doc.headings != 1: errors.append(f'{name}: expected one h1, found {doc.headings}')
        if doc.images_without_alt: errors.append(f'{name}: image without alt text')
        if doc.duplicate_ids: errors.append(f'{name}: duplicate ids: {doc.duplicate_ids}')
        for url in doc.links:
            parsed = urlsplit(url)
            if parsed.netloc and parsed.netloc != 'mathematical-tours.github.io': continue
            if parsed.scheme and parsed.scheme not in ('http', 'https'): continue
            if parsed.scheme == 'http': errors.append(f'{name}: insecure local URL: {url}')
            path = unquote(parsed.path)
            target = (root/path.lstrip('/') if path.startswith('/') else file.parent/path) if path else file
            if target.is_dir(): target /= 'index.html'
            target = target.absolute()
            checked.add(str(target))
            if not target.is_file():
                errors.append(f'{name}: missing target: {url}')
            elif parsed.fragment and target.suffix == '.html':
                # MyST builds client-side anchors, so only inspect our own static pages.
                if target.relative_to(root).as_posix() in PAGES and unquote(parsed.fragment) not in document(target).ids:
                    errors.append(f'{name}: missing anchor: {url}')
            elif target.suffix == '.pdf':
                with target.open('rb') as stream:
                    if stream.read(5) != b'%PDF-': errors.append(f'{name}: invalid PDF: {url}')
    return {'pages': len(PAGES), 'local_targets': len(checked), 'errors': errors}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('site', type=Path)
    args = parser.parse_args()
    result = check(args.site.resolve())
    print(json.dumps(result, indent=2))
    sys.exit(bool(result['errors']))

if __name__ == '__main__': main()
