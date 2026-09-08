#!/usr/bin/env python3
"""Preview the static edition at the same path used by its public build."""
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit
import argparse
import json

ROOT = Path(__file__).resolve().parents[2] / 'myst/_build/html'


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--port', type=int, default=8790)
    args = parser.parse_args()
    prefix = json.loads((ROOT / 'build-info.json').read_text())['base_url']

    class Handler(SimpleHTTPRequestHandler):
        def send_head(self):
            pathname = urlsplit(self.path).path
            if prefix and pathname in ('/', prefix):
                self.send_response(302)
                self.send_header('Location', prefix + '/')
                self.end_headers()
                return None
            if prefix and not pathname.startswith(prefix + '/'):
                self.send_error(404)
                return None
            return super().send_head()

        def translate_path(self, path):
            return super().translate_path(path[len(prefix):] if prefix else path)

    server = ThreadingHTTPServer(('127.0.0.1', args.port), partial(Handler, directory=str(ROOT)))
    print(f'Read the book at http://127.0.0.1:{args.port}{prefix}/', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()


if __name__ == '__main__':
    main()
