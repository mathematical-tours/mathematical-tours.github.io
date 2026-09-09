#!/usr/bin/env python3
"""Deploy an already-verified artifact through GitHub's Pages REST API.

Uses the runner's short-lived GITHUB_TOKEN and OIDC identity, without the
deprecated JavaScript artifact SDK bundled in actions/deploy-pages v5.
"""
import argparse
import json
import os
from pathlib import Path
import re
import signal
import time
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


class APIError(RuntimeError):
    def __init__(self, status):
        self.status = status
        super().__init__(f'GitHub request failed (HTTP {status}).')


def request_json(method, url, token, body=None):
    headers = {'Authorization': 'Bearer ' + token, 'Accept': 'application/json',
               'User-Agent': 'Mathematical-Tours-Pages'}
    if url.startswith('https://api.github.com/'):
        headers.update({'Accept': 'application/vnd.github+json', 'X-GitHub-Api-Version': '2026-03-10'})
    payload = None
    if body is not None:
        payload = json.dumps(body).encode()
        headers['Content-Type'] = 'application/json'
    try:
        with urlopen(Request(url, data=payload, headers=headers, method=method), timeout=15) as response:
            text = response.read()
            return json.loads(text) if text else {}
    except HTTPError as error:
        # Response bodies and request URLs can contain authentication material.
        raise APIError(error.code) from None
    except (URLError, TimeoutError):
        raise APIError(0) from None


def deploy(artifact_id, env, request=request_json, clock=time.monotonic,
           sleep=time.sleep, log=print, timeout=600):
    if not 0 < timeout <= 600 or artifact_id <= 0:
        raise ValueError('A positive artifact ID and a timeout of at most 600 seconds are required.')
    repository, sha = env['GITHUB_REPOSITORY'], env['GITHUB_SHA']
    if not re.fullmatch(r'[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+', repository) or not re.fullmatch(r'[0-9a-f]{40}', sha):
        raise ValueError('Invalid repository or commit identity.')
    if env.get('GITHUB_API_URL', 'https://api.github.com') != 'https://api.github.com':
        raise ValueError('This publication workflow targets GitHub.com.')
    token = env['GITHUB_TOKEN']
    identity = request('GET', env['ACTIONS_ID_TOKEN_REQUEST_URL'], env['ACTIONS_ID_TOKEN_REQUEST_TOKEN'])['value']
    endpoint = f'https://api.github.com/repos/{repository}/pages/deployments'
    # Do not retry a creation request with an uncertain outcome.
    created = request('POST', endpoint, token, {'artifact_id': artifact_id,
        'environment': 'github-pages', 'pages_build_version': sha, 'oidc_token': identity})
    deployment_id = quote(str(created.get('id') or sha), safe='')
    status_url = endpoint + '/' + deployment_id
    deadline = clock() + timeout
    finished, previous, errors = False, None, 0
    log(f'Publishing commit {sha[:12]} from verified artifact {artifact_id}.')
    try:
        while clock() < deadline:
            try:
                state = request('GET', status_url, token).get('status', 'unknown_status')
                errors = 0
            except APIError as error:
                errors += 1
                if error.status not in (0, 429, 500, 502, 503, 504) or errors >= 3:
                    raise
                log('Temporary status-request failure; checking again shortly.')
                sleep(min(15, max(0, deadline - clock())))
                continue
            if state != previous:
                log('Pages status: ' + state)
                previous = state
            if state == 'succeed':
                finished = True
                return 'https://mathematical-tours.github.io/'
            if state in ('deployment_failed', 'deployment_content_failed', 'deployment_cancelled', 'deployment_lost'):
                finished = True
                raise RuntimeError('GitHub Pages reported ' + state + '.')
            sleep(min(15, max(0, deadline - clock())))
        raise TimeoutError('GitHub Pages did not finish within its ten-minute limit.')
    finally:
        if not finished:
            try:
                request('POST', status_url + '/cancel', token)
                log('Cancelled the unfinished Pages deployment.')
            except APIError:
                log('Cancellation could not be confirmed; inspect the Pages deployment in GitHub.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--artifact-id', required=True, type=int)
    args = parser.parse_args()
    if os.environ.get('GITHUB_ACTIONS') != 'true':
        raise SystemExit('Deployment runs only in GitHub Actions with a verified workflow identity.')
    def interrupted(_number, _frame):
        raise KeyboardInterrupt('Workflow interrupted.')
    signal.signal(signal.SIGTERM, interrupted)
    try:
        url = deploy(args.artifact_id, os.environ)
    except (RuntimeError, ValueError, KeyError, TimeoutError, KeyboardInterrupt) as error:
        print('::error::' + str(error))
        raise SystemExit(1) from None
    with Path(os.environ['GITHUB_OUTPUT']).open('a') as output:
        output.write('page_url=' + url + '\n')
    print('Published successfully: ' + url)


if __name__ == '__main__':
    main()
