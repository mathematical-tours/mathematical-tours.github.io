#!/usr/bin/env python3
"""Remove reproducible build caches, retaining published PDFs and local tools."""
from pathlib import Path
import shutil

BUILD = Path(__file__).resolve().parents[1] / 'build'


def main():
    for name in ('book', 'chapters', 'figure-regeneration', 'figure-review', 'figure-processing'):
        directory = BUILD / name
        if directory.exists():
            shutil.rmtree(directory)
    (BUILD / 'build-report.json').unlink(missing_ok=True)
    print('Removed build caches. Published PDFs, Python runtime, and review evidence are retained.')


if __name__ == '__main__':
    main()
