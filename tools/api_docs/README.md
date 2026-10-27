# API reference tooling

Shared tooling that turns a repository's public API into a cross-language
reference: one page per **canonical operation** in the repo's manifest, with a
tab for each binding that exports it, showing that binding's spelling,
signatures and docstring (taken from the binding's own documentation tooling)
and any runnable examples that use the operation. Usage guides are built from
the same examples, and every example runs in CI, so the reference cannot drift
from the code.

The tooling lives only here (`tools/api_docs/`, plus the composite action
`.github/actions/api-extract`). A repository that wants an API reference keeps
**configuration and inputs only** and runs the tools from a checkout of this
repo (in CI the actions do that; see `templates/docs.yml`).

```
<repo>/docs/api/config.toml            project configuration (§0)
<repo>/api-surface.json                manifest: canonical names, tiers, per-binding spellings (§0.2)
<repo>/docs/api/groups.toml            topic grouping of the canonical operations
<repo>/docs/api/examples/<binding>/    runnable example programs, with tagged regions (§3)
<repo>/docs/api/examples/fixtures/     shared example inputs
<repo>/docs/api/guides/                usage-guide prose that pulls in example regions (§4)

tools/api_docs/<binding>.sh            extract docstrings + run that binding's examples (§1)
  -> <repo>/build/api-docs/<binding>.json   normalized docstrings (gitignored)
tools/api_docs/render.py               joins the above -> <repo>/docs/content/generated/api/ (§5)
tools/api_docs/plan.py                 prints the configured bindings as a JSON array (CI matrix)
tools/api_docs/config.py               config loader; also a CLI for the shell scripts
tools/api_docs/regions.py              example-region checker
```

Generated pages link to each other by project-root-relative path
(`/generated/api/ref/<name>`); the site assembler prefixes the project slug.

## 0. Configuration

### 0.1 `docs/api/config.toml`

```toml
manifest = "api-surface.json"      # defaults shown; all four are optional
examples = "docs/api/examples"
groups   = "docs/api/groups.toml"
guides   = "docs/api/guides"
out      = "docs/content/generated/api"

[reference]                        # all optional
title    = "EarthSciAST"           # display name (default: `title` in docs/site.toml)
summary  = "read, check and simulate the same `.esm` documents."
                                   # completes "<title> has N bindings — Julia, … — that <summary>"
repo_url = "https://github.com/EarthSciML/EarthSciAST"   # default: https://github.com/<repo in docs/site.toml>
spec           = "API_SPEC.md"     # a spec document at the repo root; with these two anchors the
spec_tiers     = "3-surface-tiers" # landing and extension pages link "(see API_SPEC.md §3)"
spec_extension = "7-extension-seams"

[bindings.<name>]                  # name: julia | typescript | python | rust | go; order = tab order
path = "..."                       # required; package/crate/module directory, relative to the repo root
label = "..."                      # display name override
examples = "..."                   # this binding's example directory (default <examples>/<name>)
display = "..."                    # package name shown on the landing page, when the manifest has no
                                   # binding_profiles.<name>.package
```

The set of `[bindings.*]` tables **is** the binding list: only those bindings
are extracted, rendered and checked. Per-binding keys:

| binding | key | meaning |
|---|---|---|
| julia | `package` | module name (`using <package>`) |
| | `env` | Julia environment with the package, `JSON3` and the examples' deps (default `<examples>/julia`) |
| | `load` | extra packages loaded before extraction (extension triggers), list |
| typescript | `package` | npm name the examples import |
| | `entry` | entry module relative to `path` (default `src/index.ts`) |
| | `dist` | built entry the examples import (default `dist/esm/index.js`) |
| python | `package` | import name |
| | `extras` | pip extras installed for the examples, list |
| | `src` | source dir relative to `path` (default `src`) |
| rust | `crate` | package name |
| | `features` | features enabled for rustdoc, list |
| | `docs_rs`, `docs_rs_features` | docs.rs root URL (default `https://docs.rs/<crate>/latest/`) and the features it builds with |
| | `nightly` | pinned nightly for rustdoc JSON (default `nightly-2026-10-01`) |
| go | `package` | import path of the package documented (`path` is its directory) |
| | `module_file` | go.mod used by CI's setup-go (default `<path>/go.mod`) |

### 0.2 The manifest (`api-surface.json`)

```json
{
  "symbols": [
    {"name": "load_path", "kind": "function", "tier": "stable",
     "bindings": {"python": ["load_path", "load_path_alias"], "rust": "load_path"},
     "notes": "optional text shown on the page"}
  ],
  "capability_profiles": {"core": {"summary": "...", "bindings": ["julia", "python"]}},
  "binding_profiles": {"python": {"package": "earthsci-ast"}},
  "planned": []
}
```

Only `symbols` is required. `kind` is `function`, `type`, `error` or
`constant`; `tier` is `stable` (rendered as reference pages) or `extension`
(listed per binding); a binding's spelling is a string or a list for an alias
set. `capability_profiles` and `binding_profiles` are optional: without
`capability_profiles` every binding is in one implicit profile and the profile
tables and notes are omitted. `binding_profiles.rust.public_modules` lists
public modules whose interiors are extension seams. `planned` lists known
harmonization gaps shown on the affected pages.

## 1. The per-binding entry point

`tools/api_docs/<binding>.sh`, for `<binding>` in `julia`, `typescript`,
`python`, `rust`, `go`, does two things for the repository named by
`--repo-root` (default: the current directory):

1. **Extracts** documentation for every spelling that the manifest lists
   under that binding, and writes `<repo>/build/api-docs/<binding>.json` (§2).
2. **Runs** every example program in that binding's example directory and
   exits non-zero if any of them fails.

`--extract-only` and `--examples-only` run one half. Each script works in CI
from a clean checkout once the language toolchain is installed, and locally
too. Python 3.11+ is needed for the tooling itself (`tomllib`); set `PYTHON`
to choose it. Build output goes under `build/` (gitignored), never `/tmp`.

In CI, `.github/actions/api-extract` does the toolchain setup (derived from the
config), the region check, the script, and uploads the artifact
`api-docs-<binding>`; `plan.py` supplies the matrix.

## 2. The normalized docstring file

```json
{
  "binding": "python",
  "package_version": "0.3.1",
  "extractor": "griffe 1.5.1",
  "symbols": [
    {
      "spelling": "load_path",
      "kind": "function",
      "signatures": ["load_path(path: str | os.PathLike) -> EsmFile"],
      "doc": "Read and parse an `.esm` file.\n\n**Parameters**\n\n- `path` — ...",
      "source": { "path": "pkg/earthsci-ast-py/src/earthsci_ast/io.py", "line": 42 },
      "native_url": null
    }
  ]
}
```

- **One entry per spelling** in the manifest for this binding, including
  both spellings where the manifest lists an alias pair. A spelling the extractor
  cannot find is still listed, with `"doc": null`, `"signatures": []` and
  `"missing": true`. Extra symbols not in the manifest are left out.
- **`kind`** uses the manifest's vocabulary: `function`, `type`, `error`,
  `constant`.
- **`signatures`** are one string per method or overload, in the binding's own
  syntax. Julia lists every method of an exported function, and Go lists a
  type's exported methods under the type's entry (as further signatures, with
  the receiver). The same goes for TypeScript/Python class members that are
  part of the public surface.
- **`doc`** is CommonMark, with the binding's own conventions converted:
  - Section syntax (Google-style `Args:`/`Returns:`, TSDoc `@param`/`@returns`,
    Rust `# Errors`, Julia `# Arguments`) becomes a bold heading followed by a
    bullet list or paragraph.
  - Cross-references (Julia `[`x`](@ref)`, Rust intra-doc ``[`Foo`]``, TSDoc
    `{@link Foo}`, Go `[Foo]`) become `[`Foo`](api:Foo)`, where `Foo` is this
    binding's spelling. The renderer resolves `api:` links to the canonical
    page, or falls back to plain code.
  - Doctest fences (`jldoctest`, Python `>>>` blocks, Rust ```` ```rust,ignore ````)
    become plain fenced code in the binding's language.
  - When a function has several docstrings (Julia methods), they are joined in
    method order, separated by a horizontal rule.
- **`source.path`** is relative to the repository root.
- **`native_url`** is a deep link to the binding's hosted docs where one exists:
  docs.rs for Rust (`https://docs.rs/earthsci-ast/latest/earthsci_ast/fn.load_path.html`),
  pkg.go.dev for Go. Otherwise `null`.

## 3. Example programs

Each binding has a directory ``<examples>/<binding>/`` of complete,
runnable programs written in that language's normal style, built against the
in-repo package (not a published release). They use the shared documents in
`<examples>/fixtures/`, located relative to the example file rather than the working
directory.

**Examples check their results.** Each step asserts what it just did, using the
language's ordinary means (`@assert`/`@test`, `assert`, `assert!`, Go
`// Output:` or `t.Fatal`, `node:assert`). An example that only runs is not
enough.

### Regions

The renderer pulls code out of the programs by region markers written as line
comments, using the language's own comment token:

```python
# region: load ops=load_path,load_string
doc = earthsci_ast.load_path(FIXTURES / "decay.esm")
assert list(doc.models) == ["Decay"]
# endregion
```

- The region id comes first. `ops=` lists the **canonical** names
  (the manifest `name`, never the binding spelling) of the operations the
  region demonstrates, each one the binding has a spelling for.
  `tools/api_docs/regions.py` is the one checker of these rules; every binding
  script and the renderer run it, and it fails on an unknown canonical name.
- Code inside the region is shown dedented. A line ending in a `hide` comment
  (`# hide`, `// hide`) runs but is not shown; use it for setup the reader does
  not need, never to hide the assertion that proves a result.
- Regions may not nest. The same id may not appear twice within one binding.

### The shared scenarios

*Region ids are a per-project convention; the table below is EarthSciAST's, shown as an example. A repo defines its own scenarios and refers to them from its guides.*

So that the guides can show one step side by side in every language, each
binding implements these region ids, for every scenario its capability profile
covers (the manifest `capability_profiles`). A binding without a stable
operation for a step leaves that region out and says why in the report. It does
not fake one.

| Region id | What it shows | Fixture |
|---|---|---|
| `load` | read a document from a path, and from a JSON string | `decay.esm` |
| `validate` | validate a good document; then show how a bad one's errors surface (result object or exception — whichever the binding does) and read the message | `decay.esm`, `invalid.esm` |
| `inspect` | list models and variables; classify them (states, algebraic unknowns, parameters) | `decay.esm` |
| `display` | render an equation or expression for people (Unicode / LaTeX / text) | `decay.esm` |
| `expressions` | build or parse an expression; free variables; substitute | — |
| `edit` | add a variable and an equation, serialize to JSON, write to a file | `decay.esm` |
| `flatten` | flatten a coupled document; list the resulting variables and equations | `coupled.esm` |
| `reactions` | derive the ODEs of a reaction system | `ozone.esm` |
| `simulate` | build a problem and solve it to t = 10 s; assert N(ln 2 / λ) ≈ 50; read the `half_life` observed (simulating bindings only) | `decay.esm` |

A binding may add further regions, named after the operation they demonstrate,
for important stable operations the scenarios do not reach.

## 4. Usage guides

A guide is a Markdown file in `docs/api/guides/` with Hugo front matter
(`title`, `description`, and `weight` to order the guides; lower comes first).
Its prose is written once for every binding. Where a step needs code, the guide
names an example region on a line of its own:

```markdown
<!-- example: load -->
```

`render.py` replaces that line with a language-tabbed block holding the region
of that id from every binding that has one, each linked to its source. A
binding without the region is left out of the block, unless the region belongs
to a capability profile the binding is outside of (for example `simulate` in
TypeScript or Go); that binding gets a one-line note saying so. The profile is
taken from the topic group (`groups.toml`) of the operations the region's
`ops=` names, and can be given explicitly with `profile=`:

```markdown
<!-- example: simulate profile=simulation -->
```

The placeholder is an HTML comment, so a guide read unrendered on GitHub shows
nothing in its place. Placeholders inside code fences are left alone. A
placeholder whose id no binding has renders a visible "not available in this
build" note and a warning.

Write guide prose that holds in every binding: describe what an operation does
and what the example shows, not one binding's types or error style. Link to
reference pages by site path, e.g. `[load_path](/generated/api/ref/load_path)`.

## 5. Rendering

```bash
python3 <earthsci.dev>/tools/api_docs/render.py --repo-root .   # Python 3.11+, standard library only
```

`render.py` reads every `build/api-docs/<binding>.json` that exists; a binding
whose file is missing still gets its pages, with a "not extracted in this
build" note in its tabs, so a local build with only some toolchains works. It
exits non-zero on a broken example region (§3), a name in `groups.toml` that is
not a stable canonical name, or a malformed docstring file, and prints warnings
for stable names in no group (shown under "Other"), placeholders with no region,
and docstring files that leave out manifest spellings. `--strict` turns the
warnings into errors. `--out`, `--api-docs`, `--examples-dir`, `--sha` and
`--repo-url` override the defaults.

Output, all under `docs/content/generated/api/` (never committed):

| Path | Page |
|---|---|
| `_index.md` | landing page: tiers, capability profiles (when defined), topics, bindings |
| `ref/<name>.md` | one per stable `(name, kind)`; `<name>-<kind>` when a name has two stable kinds |
| `topics/<group>.md` | one per `groups.toml` group |
| `bindings/<binding>.md` | a binding's stable names mapped to their canonical pages |
| `extension/<binding>.md` | a binding's extension-tier names, one anchor each |
| `coverage.md` | stable names each binding exports without a docstring |
| `guides/<guide>.md` | the guides, with placeholders filled |

`python3 regions.py --repo-root . [--binding <b>]` checks the example regions
alone, and is what a binding script calls to fail fast on a bad marker.

### Local preview

Render, then build the site from this repo's checkout with your repo's
bundle (see the site's README). The generated pages are plain Hugo markdown
with the tab markup the site's `api` layouts expect.
