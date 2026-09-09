import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

TOOLS = Path(__file__).resolve().parents[1]


def module(name):
    spec = importlib.util.spec_from_file_location(name.replace('-', '_'), TOOLS / (name + '.py'))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


prepare = module('prepare-site')
check = module('check-publication')
pages = module('deploy-pages')


class PublicationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for name, text in [('index.md', '---\ntitle: Book\n---'), ('main.pdf', '%PDF-test'),
                           ('private-data.csv', 'not published'), ('reader/index.html', '<h1>Book</h1>')]:
            path = self.root / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text)
        digest = hashlib.sha256((self.root / 'reader/index.html').read_bytes()).hexdigest()
        (self.root / 'reader/publication.json').write_text(json.dumps({'files': {'index.html': digest}}))
        self.inventory = {'templates': ['index.md'], 'assets': ['main.pdf'],
            'reader_manifest': 'reader/publication.json', 'maximum_bytes': 10000, 'maximum_files': 20}

    def test_only_declared_content_is_selected(self):
        files, assets, size = prepare.collect(self.root, self.inventory)
        self.assertEqual(files, {'index.md', 'main.pdf', 'reader/index.html', 'reader/publication.json'})
        self.assertNotIn('private-data.csv', assets)
        self.assertGreater(size, 0)

    def test_missing_download_fails_before_publication(self):
        (self.root / 'main.pdf').unlink()
        with self.assertRaisesRegex(ValueError, 'missing'):
            prepare.collect(self.root, self.inventory)

    def test_changed_reader_file_is_rejected(self):
        (self.root / 'reader/index.html').write_text('Changed outside the reader build')
        with self.assertRaisesRegex(ValueError, 'validated manifest'):
            prepare.collect(self.root, self.inventory)

    def test_symbolic_links_are_rejected(self):
        (self.root / 'main.pdf').unlink()
        (self.root / 'main.pdf').symlink_to(self.root / 'private-data.csv')
        with self.assertRaisesRegex(ValueError, 'regular file'):
            prepare.collect(self.root, self.inventory)

    def test_size_budget_is_enforced(self):
        self.inventory['maximum_bytes'] = 1
        with self.assertRaisesRegex(ValueError, 'budget'):
            prepare.collect(self.root, self.inventory)

    def test_output_audit_detects_lost_or_changed_assets(self):
        _, assets, _ = prepare.collect(self.root, self.inventory)
        (self.root / 'main.pdf').write_text('corrupt')
        (self.root / 'reader/index.html').unlink()
        result = check.check(self.root, dict(self.inventory, assets=assets))
        self.assertEqual(len(result['errors']), 2)


class DeploymentTests(unittest.TestCase):
    def setUp(self):
        self.env = {'GITHUB_REPOSITORY': 'owner/site', 'GITHUB_SHA': 'a' * 40,
            'GITHUB_TOKEN': 'secret-api-token', 'ACTIONS_ID_TOKEN_REQUEST_URL': 'https://oidc.example/token',
            'ACTIONS_ID_TOKEN_REQUEST_TOKEN': 'secret-runner-token'}
        self.calls, self.logs, self.states = [], [], ['succeed']
        self.now = 0

    def request(self, method, url, token, body=None):
        self.calls.append((method, url, token, body))
        if url == self.env['ACTIONS_ID_TOKEN_REQUEST_URL']:
            return {'value': 'secret-oidc-identity'}
        if method == 'POST':
            return {} if url.endswith('/cancel') else {'id': 'deployment-id'}
        state = self.states.pop(0) if len(self.states) > 1 else self.states[0]
        if isinstance(state, Exception):
            raise state
        return {'status': state}

    def sleep(self, seconds):
        self.now += seconds

    def deploy(self, **kwargs):
        return pages.deploy(123, self.env, request=self.request, clock=lambda: self.now,
                            sleep=self.sleep, log=self.logs.append, **kwargs)

    def test_waits_for_success_and_keeps_credentials_out_of_logs(self):
        self.states = ['queued', 'syncing_files', 'succeed']
        self.assertEqual(self.deploy(), 'https://mathematical-tours.github.io/')
        post = [c for c in self.calls if c[0] == 'POST']
        self.assertEqual(len(post), 1)
        self.assertEqual(post[0][3]['artifact_id'], 123)
        self.assertEqual(post[0][3]['oidc_token'], 'secret-oidc-identity')
        self.assertNotIn('secret-', '\n'.join(self.logs))

    def test_timeout_cancels_the_incomplete_deployment(self):
        self.states = ['syncing_files']
        with self.assertRaises(TimeoutError):
            self.deploy(timeout=30)
        self.assertTrue(self.calls[-1][1].endswith('/cancel'))

    def test_terminal_server_failure_is_not_reported_as_success(self):
        self.states = ['deployment_content_failed']
        with self.assertRaisesRegex(RuntimeError, 'deployment_content_failed'):
            self.deploy()
        self.assertFalse(self.calls[-1][1].endswith('/cancel'))

    def test_transient_read_error_recovers_without_creating_twice(self):
        self.states = [pages.APIError(503), 'succeed']
        self.deploy()
        self.assertEqual(len([c for c in self.calls if c[0] == 'POST']), 1)

    def test_authorization_failure_cancels_and_fails(self):
        self.states = [pages.APIError(403)]
        with self.assertRaises(pages.APIError):
            self.deploy()
        self.assertTrue(self.calls[-1][1].endswith('/cancel'))

    def test_unknown_creation_outcome_is_not_retried(self):
        calls = []
        def failed_create(method, *args):
            calls.append(method)
            if method == 'GET':
                return {'value': 'secret-oidc-identity'}
            raise pages.APIError(0)
        with self.assertRaises(pages.APIError):
            pages.deploy(123, self.env, request=failed_create)
        self.assertEqual(calls, ['GET', 'POST'])

    def test_timeout_cannot_exceed_github_limit(self):
        with self.assertRaises(ValueError):
            self.deploy(timeout=601)
        self.assertEqual(self.calls, [])


if __name__ == '__main__':
    unittest.main()
