# earthsci.dev design

One Hugo site (`earthsci.dev`) holds the documentation of six repos. This repo
also owns the shared API-doc tooling, which the other repos run from here
instead of carrying copies.

| id | repo | notes |
|---|---|---|
| `ast` | EarthSciML/EarthSciAST | julia, typescript, python, rust, go bindings; narrative pages |
| `io` | EarthSciML/EarthSciIO | python, julia, rust |
| `astdiff` | EarthSciML/EarthSciASTDiff | julia, python, rust |
| `integrator` | EarthSciML/EarthSciIntegrator | julia, rust |
| `discretizations` | EarthSciML/EarthSciDiscretizations | generated catalog, no API ref |
| `models` | EarthSciML/EarthSciModels | generated catalog, no API ref |

Site URLs: `earthsci.dev/` (intro, quick start), `earthsci.dev/<id>/...`.

## Data flow

```
source repo (push to main)
  docs.yml (small caller workflow)
    on push to main (version `dev`) and on push of a tag vX.Y.Z (that version)
    job extract (matrix: bindings)   uses EarthSciML/earthsci.dev/.github/actions/api-extract@main
    job bundle                       project prebuild steps, then
                                     uses EarthSciML/earthsci.dev/.github/actions/docs-bundle@main
        render API ref -> docs/content/generated/api/
        assemble bundle (site.toml + content/ + static/ + api/<version>/)
        merge bundle into orphan branch `docs-bundle` of the SOURCE repo
          (keeps every version's API tree; tools/publish.py)
        repository_dispatch(event_type=docs-updated, payload {id, repo, sha}) -> EarthSciML/earthsci.dev
earthsci.dev (.github/workflows/site.yml, on repository_dispatch / push / manual)
  tools/assemble.py: for every project in projects.toml, fetch `docs-bundle`
  branch of its repo, copy into content/<id>/ and static/<id>/, place the
  versioned API trees (see Versions), rewrite root-relative links to
  /<id>/..., write data/projects/<id>.json
  hugo build, link check, deploy to GitHub Pages (CNAME earthsci.dev)
```

Every build refetches every project's latest `docs-bundle`, so the result is
the same whichever repo triggered it; the dispatch payload only names what
changed (for the run summary). A repo with no `docs-bundle` branch yet is
skipped with a warning. On pull requests the source repo's workflow runs
extract + render strictly and builds the bundle but does not push or dispatch.

Secrets: each source repo needs `DOCS_APP_ID` and `DOCS_APP_PRIVATE_KEY` (a
GitHub App installed on EarthSciML/earthsci.dev with contents:write; org-level
secrets work). The bundle job mints a short-lived token from them with
`actions/create-github-app-token` and passes it to `docs-bundle` as
`dispatch-token`. Pushing `docs-bundle` uses the repo's own `GITHUB_TOKEN` (needs
`contents: write` in the bundle job).

## Versions

Only the API reference (`/<id>/generated/api/**`) is versioned; the narrative
docs always come from main. Versions are `dev` (main) and one `vX.Y.Z` per
release tag (`vX.Y.Z-pre` tags are versions too, but never `latest`).

```
/<id>/generated/api/            latest: the newest non-prerelease tag (dev while there is none)
/<id>/generated/api/dev/        main
/<id>/generated/api/v1.2.0/     a release
```

- **Bundle layout.** A bundle is `bundle.json`, `content/`, `static/` (main's
  narrative docs; the API reference is *not* in `content/`) plus
  `api/<version>/` (the rendered `docs/content/generated/api` tree for one
  version) and `versions.json` (`{"latest", "versions"}`). The set of versions is
  just the directories under `api/`, so there is no separate state.
- **Source repo.** `docs.yml` runs on pushes to main and on tags `v*`. The
  `docs-bundle` action derives the version from the ref (override with its
  `version` input), runs `tools/bundle.py --version <v>` (a tag bundle carries only
  `api/<tag>/`, built from the code at that tag with the docs config at that tag)
  and `tools/publish.py`, which clones the existing `docs-bundle` branch, merges
  (`dev` replaces content/static/`api/dev`; a tag adds or replaces `api/<tag>`),
  prunes, and pushes the result as a single orphan commit (so history does not
  grow) with `--force-with-lease`, retrying after a conflicting concurrent publish.
  Tags that are not `vX.Y.Z[-pre]` are ignored with a notice.
- **Bounded size.** Only the newest patch of each of the newest N minor series
  is kept (plus `dev`, and the newest prerelease if it is newer than `latest`).
  N is `max-minor-versions` on the action (default 8) on the branch, and
  `max_minor_versions` in `projects.toml` (default 8, per-project override) on the
  site. 0 keeps everything.
- **Site.** `tools/assemble.py` copies the `latest` tree to `generated/api/`
  and each kept version to `generated/api/<version>/`. Links inside a version tree
  (`/generated/api/...`) are rewritten to stay in that version, then prefixed
  with `/<id>`; the `latest` copy keeps unversioned links. `data/projects/<id>.json`
  gets `api = {latest, latest_is_release, versions}`.
- **UI.** `layouts/partials/api-version.html` (included by `baseof.html` on API
  pages) shows a version dropdown (`latest (v1.2.0)`, `dev`, each version) and a
  banner on `dev`, older and pre-release versions linking to the latest release.
  A little JS retargets the dropdown links to the same page in the other version
  when it exists. Only the default (`latest`) copy is indexed by Pagefind.
- **No release yet.** The project has `dev` only; `latest` is `dev` and the
  dropdown lists just `dev`.
- **Caveat.** Per-symbol source links in generated pages point at the commit that
  built that version, but the `[spec]` links render.py writes (`blob/main/...`)
  point at main in every version.

## Source-repo docs contract

A source repo keeps only content and configuration; no `hugo.toml`, no layouts,
no extraction scripts.

### `docs/site.toml` (required)

```toml
id = "io"                           # URL slug; must match projects.toml here
title = "EarthSciIO"
description = "One-sentence summary shown on the home page and section index."
repo = "EarthSciML/EarthSciIO"
content = "docs/content"            # default; Hugo markdown, `_index.md` is the section landing page
static = "docs/static"              # default; optional; copied to /<id>/...
prebuild = ["python3 docs/generate_catalog.py"]   # optional; informational here, the
                                    # caller workflow runs them before docs-bundle
[[nav]]                             # optional; project header links. paths are relative to /<id>/
title = "Operators"
path = "operators/"
```

The API reference link (`generated/api/`) is added to the nav automatically
when `docs/api/config.toml` exists. Content links that are root-relative
(`/generated/api/ref/x`, `/plots/a.png`) are written as if the project were the
site root; assemble prefixes `/<id>`. Hugo `title`/`weight` front matter works
as usual. Generated content (`docs/content/generated/`) is never committed.

### `docs/api/config.toml` (only for projects with an API reference)

```toml
manifest = "api-surface.json"       # default; canonical-operation manifest (see tools/api_docs/README.md)
examples = "docs/api/examples"      # default
groups   = "docs/api/groups.toml"   # default
guides   = "docs/api/guides"        # default

[reference]                         # prose for the generated landing pages
summary = "..."                     # the project's one-paragraph description of its bindings

[bindings.python]                   # binding names: julia typescript python rust go
path = "pkg/earthsci-ast-py"        # package root, relative to the repo root
package = "earthsci_ast"            # importable/crate/module name (per-binding extras documented in tools/api_docs/README.md)
```

The extractors, the example runner scripts, `regions.py` and `render.py` take
`--repo-root` and read this file; nothing project-specific is hard-coded in
`tools/api_docs/`. Build output goes under `<repo>/build/` (gitignored).

### `api-surface.json`

The manifest format is the one EarthSciAST already uses (`symbols` with
`name`, `kind`, `tier`, `bindings`; plus `capability_profiles`,
`binding_profiles`). New repos may omit `capability_profiles`/`binding_profiles`
(then every binding is in one `core` profile) and may mark everything `stable`.

## This repo

```
hugo.toml, layouts/, static/css/   the shared site theme (no external theme)
content/_index.md, content/quickstart/ ...   intro and quick start (owned here)
projects.toml                      [[project]] id, repo, title, order
tools/api_docs/                    generalized API-doc tooling + README
tools/assemble.py                  bundle fetch/copy/link-rewrite (+ tests)
tools/bundle.py                    builds a bundle from a source repo (used by docs-bundle)
tools/publish.py, tools/versions.py  merge a bundle into the docs-bundle branch; version ordering/pruning
.github/actions/api-extract        composite action: one binding's extract + examples
.github/actions/docs-bundle        composite action: render + bundle + push + dispatch
.github/workflows/site.yml         build and deploy
templates/docs.yml                 the caller workflow a source repo installs
```

## Cross-project links

`assemble.py` prefixes every root-relative link with `/<id>`, so a link from one
project to another must be an absolute URL (`https://earthsci.dev/ast/...`).
