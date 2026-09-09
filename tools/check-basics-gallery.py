#!/usr/bin/env python3
"""Check the bilingual figure gallery against the books and illustration inputs."""
import hashlib
import json
from pathlib import Path


def check(root):
    data=json.loads((root/'_data/basics_figures.json').read_text())
    errors=[]
    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None
    for name,sha in data['inputs'].items():
        if digest(root/'book-basics-sources'/name)!=sha:
            errors.append('Illustration input changed; rebuild the gallery: '+name)
    for lang,book in data['books'].items():
        if digest(root/book['pdf'].lstrip('/'))!=book['sha256']:
            errors.append('The '+lang+' book changed; refresh figure page links.')
    for lang in ('en','fr'):
        rows=[f for f in data['figures'] if f['lang']==lang]
        numbers=[tuple(map(int,f['number'].split('.'))) for f in rows]
        if len(rows)!=data['count'] or numbers!=sorted(set(numbers)):
            errors.append('Missing, duplicate, or unordered figures in '+lang)
        for row in rows:
            for key in ('pdf','jpg'):
                path=root/row[key].lstrip('/')
                if digest(path)!=row[key+'_sha256']:
                    errors.append('Missing or modified download: '+row[key])
            if not (root/row['thumbnail'].lstrip('/')).is_file():errors.append('Missing preview: '+row['id'])
            if row['book']!=data['books'][lang]['pdf']:errors.append('Wrong edition link: '+row['id'])
    return {'figures_per_edition':data['count'],'editions':len(data['books']),
            'pdf_downloads':len(data['figures']),'jpg_downloads':len(data['figures']),'errors':errors}


if __name__=='__main__':
    result=check(Path(__file__).resolve().parents[1])
    print(json.dumps(result,indent=2))
    raise SystemExit(bool(result['errors']))
