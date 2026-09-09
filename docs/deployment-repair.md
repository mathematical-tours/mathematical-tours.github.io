# Pages deployment repair — 9 September 2026

The artwork update in commit `b6ed18fd` built successfully in GitHub Actions but
[its deployment](https://github.com/mathematical-tours/mathematical-tours.github.io/actions/runs/34351731751)
remained in `syncing_files` until the ten-minute timeout. The upload archive was
535,386,173 bytes. The logs do not establish why the Pages service stopped making
progress. This was a deployment failure, not a failed Jekyll build.

The generated workflow also used the Node 20 upload action and emitted DEP0040
(`punycode`) and DEP0169 (`url.parse()`) warnings. The deployment action’s deprecated
artifact dependency is a [known upstream issue](https://github.com/actions/deploy-pages/issues/413).
These warnings were separate from the synchronization timeout.

## Changes

- Added a versioned publishing workflow and an explicit publication inventory.
  The prepared source is about 489 MB. It contains the current website, all linked
  PDFs/notebooks, compatibility downloads, and the complete validated interactive
  book. Unused datasets, book source trees, intermediate TeX files and unlinked
  auxiliary materials remain available in Git without being uploaded to Pages.
- Before upload, check the pages and local links, every asset checksum, file types,
  symbolic/hard links, total file count, and size. The budget is 600 MB and 5,000
  files, below GitHub’s 1 GB publication limit. The full interactive edition is
  preserved byte for byte; no figure or PDF quality is reduced.
- Pin official checkout and upload actions to Node 24 releases. Upload the regular
  Pages tar artifact with the current uploader, avoiding the outdated uploader in
  the generated branch-based workflow.
- Deploy through GitHub’s documented Pages REST API using the standard short-lived
  workflow token and OIDC identity. No JavaScript dependency or warning suppression
  is needed for this step. Success requires the server’s `succeed` status; transient
  read errors are retried, creation requests are never blindly repeated, and
  unfinished deployments are cancelled at the ten-minute limit or on interruption.
- Add regression checks for missing/altered downloads, unexpected publication
  inputs, symbolic links, size limits, credential-safe logging, transient API
  failures, terminal failures, cancellation and uncertain creation outcomes.

## Reproduction

From the website repository root:

```sh
python3 -m unittest discover -s tools/tests -v
python3 tools/prepare-site.py --destination .pages-source --report .pages-inputs.json
bundle exec jekyll build --source .pages-source --destination _site
python3 tools/check-site.py _site
python3 tools/check-publication.py _site --source-map .pages-inputs.json
```

GitHub Pages uses **GitHub Actions** as its build source; the public URL and
`master` branch stay the same. The workflow’s concurrency group serializes
publication so another push cannot interrupt a deployment halfway through.
