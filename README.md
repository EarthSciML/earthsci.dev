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
```

Python 3.11+ is needed (`tomllib`); on older versions `pip install tomli`.

## Secrets

- `SITE_DISPATCH_TOKEN` (in each source repo, or org-wide): token with
  `contents: write` on this repo, used to send the `docs-updated` dispatch.
- `DOCS_FETCH_TOKEN` (optional, here): read access to the source repos if any
  is private. Public repos need nothing.

Pages must be set to deploy from GitHub Actions, with the custom domain in
`CNAME`.
