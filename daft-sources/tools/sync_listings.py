#!/usr/bin/env python3
"""Derive TeX listing slices from the actual, tested Python functions."""
import ast
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]

def main():
    manifest=json.loads((ROOT/'docs/listing-manifest.json').read_text())
    folder=ROOT/'listings/tex'
    folder.mkdir(parents=True,exist_ok=True)
    count=0
    for item in manifest:
        if item['status']!='python':
            continue
        path=ROOT/'listings'/item['python_file']
        tree=ast.parse(path.read_text())
        nodes={n.name:n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
        slices=[]
        for name in item['functions']:
            n=nodes[name]
            slices.append(r'\PythonListing{'+item['python_file']+'}{'+str(n.lineno)+'}{'+str(n.end_lineno)+'}')
        (folder/(item['id']+'.tex')).write_text('\n'.join(slices)+'\n')
        count+=1
    print(f'Synchronized {count} Python listings.')

if __name__=='__main__':
    main()
