#!/usr/bin/env python3
"""Run the mathematical checks and verify execution of every printed function."""
import ast
import json
from pathlib import Path
import platform
import sys
import time
import unittest
from importlib.metadata import version

ROOT = Path(__file__).resolve().parents[1]


def main():
    sys.path.insert(0, str(ROOT / 'listings'))
    manifest = json.loads((ROOT / 'docs/listing-manifest.json').read_text())
    called = set()
    prefix = str(ROOT / 'listings') + '/'

    def record(frame, event, arg):
        if event == 'call' and frame.f_code.co_filename.startswith(prefix):
            called.add((Path(frame.f_code.co_filename).name, frame.f_code.co_name))

    suite = unittest.defaultTestLoader.discover(str(ROOT / 'tests'))
    start = time.monotonic()
    sys.setprofile(record)
    try:
        result = unittest.TextTestRunner(verbosity=2).run(suite)
    finally:
        sys.setprofile(None)
    report = []
    for item in manifest:
        filename = item['python_file']
        tree = ast.parse((ROOT / 'listings' / filename).read_text())
        functions = {n.name for n in tree.body if isinstance(n, ast.FunctionDef)}
        missing = [name for name in item['functions']
                   if name not in functions or (filename, name) not in called]
        report.append({'id': item['id'], 'module': filename,
                       'functions': item['functions'], 'unexecuted': missing})
    missing = [item for item in report if item['unexecuted']]
    evidence = {
        'python': platform.python_version(),
        'dependencies': {name: version(name) for name in
                         ['numpy', 'scipy', 'sympy', 'galois', 'matplotlib', 'Pillow']},
        'tests_run': result.testsRun, 'failures': len(result.failures),
        'errors': len(result.errors), 'elapsed_seconds': round(time.monotonic()-start, 2),
        'listings': len(report), 'listings_fully_executed': len(report)-len(missing),
        'scope': 'Every displayed function was called by the independent mathematical checks; this is not a claim of exhaustive branch coverage.',
        'programs': report,
    }
    (ROOT / 'docs/program-verification.json').write_text(json.dumps(evidence, indent=2)+'\n')
    print(f"{evidence['listings_fully_executed']}/{len(report)} printed programs fully exercised.")
    for item in missing:
        print(item['id'] + ': ' + ', '.join(item['unexecuted']))
    return 0 if result.wasSuccessful() and not missing else 1


if __name__ == '__main__':
    raise SystemExit(main())
