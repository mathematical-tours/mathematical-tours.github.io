#!/usr/bin/env python3
"""Verify the unchanged pre-revision archive against its original SHA-256 manifest."""
from pathlib import Path
import hashlib
import json
import unicodedata

ROOT = Path(__file__).resolve().parents[1]
manifest = json.loads((ROOT/'legacy/original-manifest.json').read_text())
archive = ROOT/'legacy'/manifest['snapshot']
# Git stores precomposed names on macOS; compare canonically equivalent paths.
nfc = lambda name: unicodedata.normalize('NFC', name)
expected = {nfc(record['path']):record for record in manifest['files']}
actual = {nfc(str(p.relative_to(archive))):p for p in archive.rglob('*') if p.is_file()}
assert actual.keys() == expected.keys(), ('Archive file inventory differs', actual.keys() ^ expected.keys())
for name, record in expected.items():
    data = actual[name].read_bytes()
    assert len(data) == record['bytes'], name
    assert hashlib.sha256(data).hexdigest() == record['sha256'], name
print(f'{len(expected)} original files verified byte for byte.')
