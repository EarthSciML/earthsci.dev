# earthsci.dev

Source of [earthsci.dev](https://earthsci.dev), the documentation site for the
EarthSciML projects (EarthSciAST, EarthSciIO, EarthSciASTDiff,
EarthSciIntegrator, EarthSciDiscretizations, EarthSciModels). It also owns the
shared API-documentation tooling the other repos run from here.

Each project keeps its docs (content, config) in its own repo. On every push to
`main`, the repo publishes a docs bundle and notifies this repo, which
re-assembles all bundles and redeploys the whole site. See
[DESIGN.md](DESIGN.md) for the architecture and the contract a source repo
follows.

## Layout

| Path | Purpose |
|---|---|
| `content/` | The home page and quick start (owned here). Project content is assembled into `content/<id>/` at build time and is not committed. |
| `hugo.toml`, `layouts/`, `static/css/` | The shared theme |
| `projects.toml` | The projects in the site |
| `tools/assemble.py`, `tools/bundle.py` | Bundle fetching/assembly and bundle building (tests in `tests/`) |
| `tools/publish.py`, `tools/versions.py` | Merging a bundle into a repo's `docs-bundle` branch; API version ordering and pruning (`DESIGN.md`, "Versions") |
| `tools/api_docs/` | Generalized API-reference extraction and rendering |
| `.github/actions/` | Composite actions source repos use (`api-extract`, `docs-bundle`) |
| `.github/workflows/site.yml` | Site build and deploy |
| `templates/docs.yml` | The workflow a source repo installs |

## Building locally

```bash
python3 -m unittest discover -s tests          # tool tests
python3 tools/assemble.py                      # fetch every project's docs-bundle branch
python3 tools/assemble.py --local io=/path/to/bundle --only io   # or use a local bundle
hugo server                                    # Hugo extended 0.131 or newer
```

To produce a bundle from a source checkout (after any generation steps listed
in its `docs/site.toml`):

```bash
python3 tools/bundle.py --repo-root /path/to/repo --out build/bundles/<id>
python3 tools/bundle.py --repo-root /path/to/repo --out build/bundles/<id>-v1 --version v1.0.0   # a release's API tree only
python3 tools/publish.py --bundle build/bundles/<id>-v1 --version v1.0.0 --existing build/bundles/<id>-merged   # merge locally (no push)
```

The API reference is versioned (`dev`, `latest`, `vX.Y.Z`); see `DESIGN.md`,
"Versions". To preview several versions, merge a `dev` bundle and one or more
tag bundles into one directory with `publish.py --existing`, and pass that
directory to `assemble.py --local`.

Python 3.11+ is needed (`tomllib`); on older versions `pip install tomli`.

## Secrets

- `DOCS_APP_ID` and `DOCS_APP_PRIVATE_KEY` (in each source repo, or org-wide):
  a GitHub App installed on this repo with `contents: write`; the source
  repos' workflows mint a token from it to send the `docs-updated` dispatch.
- `DOCS_FETCH_TOKEN` (optional, here): read access to the source repos if any
  is private. Public repos need nothing.

Pages must be set to deploy from GitHub Actions, with the custom domain in
`CNAME`.
