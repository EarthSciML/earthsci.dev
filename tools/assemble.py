#!/usr/bin/env python3
"""Assemble every project's docs bundle into this Hugo site.

For each [[project]] in projects.toml: obtain its bundle (the `docs-bundle`
branch of its repo, or a local directory given with --local), copy it into
content/<id>/ and static/<id>/, prefix root-relative links with /<id>, and write
data/projects/<id>.json for the header. A project with no bundle is skipped with
a warning.

Docs are versioned as a whole: the bundle holds one complete tree per version
(`dev`, `vX.Y.Z`; see versions.py). The newest release (or `dev` when there is
none) is `latest`, served at /<id>/; every kept version is served again at
/<id>/<version>/, with its links pointing at that same version. Only `latest`
contributes to the site-wide taxonomies, so the versioned copies lose their
taxonomy (and url/alias) front matter. A bundle in the older flat layout, or a
single tree built by bundle.py, is treated as one version.

    python3 tools/assemble.py
    python3 tools/assemble.py --local discretizations=/path/to/bundle

Set GH_TOKEN to fetch bundles from private repos.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

try:
    import tomllib
except ImportError:  # Python < 3.11
    import tomli as tomllib  # type: ignore

sys.path.insert(0, str(Path(__file__).resolve().parent))
import versions  # noqa: E402

BRANCH = "docs-bundle"
DEFAULT_MAX_MINOR = 5
# Front matter that must not repeat in a versioned copy besides the taxonomies.
VERSIONED_DROP = ("url", "aliases")
TEXT_SUFFIXES = {".md", ".html"}

_FENCE = re.compile(r"^\s{0,3}(`{3,}|~{3,})")
_CODE_SPAN = re.compile(r"(`+[^`\n]*`+)")
# [text](/path   ![alt](/path   [ref]: /path   href="/path"   src='/path'
_LINK_PATTERNS = [
    re.compile(r"(\]\()(/(?!/))"),
    re.compile(r"(^\s{0,3}\[[^\]\n]+\]:\s*)(/(?!/))"),
    re.compile(r"""(\b(?:href|src|poster|action|data-src)\s*=\s*["'])(/(?!/))""", re.I),
    # {{< ref "/x" >}}  {{< relref "/x" >}}  (also {{% ... %}})
    re.compile(r"""(\{\{[<%]\s*(?:rel)?ref\s+["'])(/(?!/))"""),
]


def _prefixer(project_id: str, skip: set[str], base: str):
    def repl(m: re.Match, text: str) -> str:
        rest = text[m.end():]
        first = re.split(r"[/?#)\"'\s]", rest, maxsplit=1)[0]
        if first == project_id or first in skip:
            return m.group(0)
        return f"{m.group(1)}{base}/"
    return repl


def rewrite_links(text: str, project_id: str, skip: set[str] = frozenset(),
                  version: str | None = None) -> str:
    """Prefix root-relative links with /<project_id> (/<project_id>/<version> for
    a versioned copy), leaving code alone."""
    base = f"/{project_id}/{version}" if version else f"/{project_id}"
    repl = _prefixer(project_id, set(skip), base)
    patterns = [(p, repl) for p in _LINK_PATTERNS]
    out: list[str] = []
    fence: str | None = None
    for line in text.splitlines(keepends=True):
        m = _FENCE.match(line)
        if fence is None and m:
            fence = m.group(1)[0] * 3
            out.append(line)
            continue
        if fence is not None:
            if m and m.group(1).startswith(fence):
                fence = None
            out.append(line)
            continue
        parts = _CODE_SPAN.split(line)
        for i in range(0, len(parts), 2):  # even indexes are outside code spans
            seg = parts[i]
            for pat, fn in patterns:
                seg = pat.sub(lambda mm, seg=seg, fn=fn: fn(mm, seg), seg)
            parts[i] = seg
        out.append("".join(parts))
    return "".join(out)


def rewrite_tree(root: Path, project_id: str, skip: set[str],
                 version: str | None = None) -> int:
    n = 0
    for p in root.rglob("*"):
        if p.suffix in TEXT_SUFFIXES and p.is_file():
            old = p.read_text(encoding="utf-8")
            new = rewrite_links(old, project_id, skip, version)
            if new != old:
                p.write_text(new, encoding="utf-8")
                n += 1
    return n


def fetch_bundle(repo: str, dest: Path) -> bool:
    """Shallow-clone the docs-bundle branch of `repo`; False if it has none."""
    token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
    auth = f"x-access-token:{token}@" if token else ""
    url = f"https://{auth}github.com/{repo}.git"
    if dest.exists():
        shutil.rmtree(dest)
    r = subprocess.run(
        ["git", "clone", "--quiet", "--depth", "1", "--branch", BRANCH, url, str(dest)],
        capture_output=True, text=True,
        env={**os.environ, "GIT_TERMINAL_PROMPT": "0"},
    )
    if r.returncode != 0:
        print(f"assemble.py: warning: no `{BRANCH}` branch for {repo}: "
              f"{r.stderr.strip().splitlines()[-1] if r.stderr.strip() else 'clone failed'}",
              file=sys.stderr)
        return False
    return True


def _clean(site: Path, pid: str) -> None:
    for sub in (f"content/{pid}", f"static/{pid}", f"data/projects/{pid}.json"):
        p = site / sub
        if p.is_dir():
            shutil.rmtree(p)
        elif p.exists():
            p.unlink()


def strip_front_matter_keys(text: str, keys) -> str:
    """Drop the top-level `keys` from a page's YAML (---) or TOML (+++) front
    matter, with their (indented or `- `) continuation lines."""
    keys = set(keys)
    m = re.match(r"(---|\+\+\+)[ \t]*\r?\n", text)
    if not m or not keys:
        return text
    delim = m.group(1)
    end = re.search(r"^%s[ \t]*\r?$" % re.escape(delim), text[m.end():], re.M)
    if not end:
        return text
    head = text[m.end():m.end() + end.start()]
    rest = text[m.end() + end.start():]
    out: list[str] = []
    skipping = False
    depth = 0
    in_tables = False
    for line in head.splitlines(keepends=True):
        if delim == "---":
            top = re.match(r"([A-Za-z_][\w-]*)\s*:", line)
            if top:
                skipping = top.group(1) in keys
            elif skipping and not (line[:1] in " \t" or line.startswith("- ")
                                   or line.strip() in ("", "-")):
                skipping = False
        else:
            if depth > 0:
                depth += line.count("[") - line.count("]")
                continue
            if line.lstrip().startswith("["):
                in_tables = True
            top = None if in_tables else re.match(r"\s*([A-Za-z_][\w-]*)\s*=", line)
            skipping = bool(top and top.group(1) in keys)
            if skipping:
                depth = max(0, line.count("[") - line.count("]"))
        if not skipping:
            out.append(line)
    return text[:m.end()] + "".join(out) + rest


def strip_tree(root: Path, keys) -> int:
    n = 0
    for p in root.rglob("*.md"):
        old = p.read_text(encoding="utf-8")
        new = strip_front_matter_keys(old, keys)
        if new != old:
            p.write_text(new, encoding="utf-8")
            n += 1
    return n


def taxonomy_keys(site: Path) -> set[str]:
    """Front matter keys that feed the site-wide Hugo taxonomies (singular and
    plural names), from hugo.toml."""
    cfg_path = site / "hugo.toml"
    if not cfg_path.is_file():
        return set()
    with cfg_path.open("rb") as f:
        tax = tomllib.load(f).get("taxonomies", {})
    return set(tax) | set(tax.values())


def version_trees(bundle: Path, scratch: Path) -> dict[str, Path]:
    """The version trees in a bundle: {version: tree directory}.

    A bundle is a `docs-bundle` checkout (one directory per version), a bundle in
    the older flat layout (main plus api/<version>/), or a single tree from
    tools/bundle.py, which is the version its bundle.json names (default `dev`)."""
    found = versions.scan(bundle)
    if found:
        return {n: bundle / n for n in found}
    if not (bundle / "bundle.json").is_file():
        return {}
    flat = scratch / "flat"
    shutil.copytree(bundle, flat, ignore=shutil.ignore_patterns("api"))
    if (bundle / "api" / versions.DEV).is_dir():
        (flat / "content" / "generated").mkdir(parents=True, exist_ok=True)
        shutil.copytree(bundle / "api" / versions.DEV, flat / "content" / "generated" / "api",
                        dirs_exist_ok=True)
    name = json.loads((bundle / "bundle.json").read_text(encoding="utf-8")).get("version", "")
    return {name if versions.is_version(name) else versions.DEV: flat}


def has_bundle(bundle: Path) -> bool:
    return bool(versions.scan(bundle)) or (bundle / "bundle.json").is_file()


def place_tree(site: Path, project: dict, tree: Path, skip: set[str], drop: set[str],
               version: str | None) -> dict:
    """Copy one version tree into content/<id>[/<version>] and static/<id>[/<version>]
    with links rewritten; returns its bundle.json. `version` None is `latest`."""
    pid = project["id"]
    meta = json.loads((tree / "bundle.json").read_text(encoding="utf-8"))
    if meta["id"] != pid:
        raise SystemExit(f"assemble.py: bundle for {project['repo']} has id "
                         f"{meta['id']!r}, projects.toml says {pid!r}")
    sub = f"{pid}/{version}" if version else pid
    content = site / "content" / sub
    content.mkdir(parents=True, exist_ok=True)
    if (tree / "content").is_dir():
        shutil.copytree(tree / "content", content, dirs_exist_ok=True)
    if (tree / "static").is_dir():
        shutil.copytree(tree / "static", site / "static" / sub, dirs_exist_ok=True)
    if version:
        strip_tree(content, drop)
    rewrite_tree(content, pid, skip, version)
    index = content / "_index.md"
    if not index.exists():
        index.write_text(
            f'---\ntitle: "{meta["title"]}"\ndescription: "{meta["description"]}"\n---\n',
            encoding="utf-8",
        )
    return meta


def assemble_project(site: Path, project: dict, bundle: Path, skip: set[str],
                     max_minor: int = DEFAULT_MAX_MINOR, drop: set[str] = frozenset()) -> dict:
    pid = project["id"]
    _clean(site, pid)
    with tempfile.TemporaryDirectory(prefix="assemble-") as scratch:
        trees = version_trees(bundle, Path(scratch))
        kept = versions.order(versions.prune(list(trees), project.get("max_minor_versions", max_minor)))
        top = versions.effective_latest(kept)
        # `latest` first: the versions go in subdirectories of its tree.
        top_meta = place_tree(site, project, trees[top], skip, set(), None)
        metas = {name: place_tree(site, project, trees[name], skip, set(drop) | set(VERSIONED_DROP), name)
                 for name in kept}
    data = {
        "id": pid,
        "title": top_meta.get("title") or project["title"],
        "description": top_meta["description"],
        "repo": top_meta.get("repo") or project["repo"],
        "order": project.get("order", 100),
        "nav": top_meta.get("nav", []),
        "sha": top_meta.get("sha", ""),
        "versions": {
            "latest": top,
            "latest_is_release": bool(top and versions.is_release(top)),
            "list": kept,
        },
        "nav_by_version": {name: m.get("nav", []) for name, m in metas.items()},
    }
    (site / "data" / "projects").mkdir(parents=True, exist_ok=True)
    (site / "data" / "projects" / f"{pid}.json").write_text(
        json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return data


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--site-root", type=Path, default=Path(__file__).resolve().parent.parent)
    ap.add_argument("--local", action="append", default=[], metavar="ID=PATH",
                    help="use the bundle directory PATH for project ID instead of fetching")
    ap.add_argument("--only", action="append", default=[], metavar="ID",
                    help="assemble only these projects")
    ap.add_argument("--workdir", type=Path, help="where fetched bundles go (default: <site>/build/bundles)")
    args = ap.parse_args(argv)

    site = args.site_root
    with (site / "projects.toml").open("rb") as f:
        cfg = tomllib.load(f)
    skip = set(cfg.get("global_prefixes", []))
    drop = taxonomy_keys(site)
    max_minor = cfg.get("max_minor_versions", DEFAULT_MAX_MINOR)
    local = dict(item.split("=", 1) for item in args.local)
    workdir = args.workdir or site / "build" / "bundles"

    done, skipped = [], []
    for project in cfg["project"]:
        pid = project["id"]
        if args.only and pid not in args.only:
            continue
        if pid in local:
            bundle = Path(local[pid])
            ok = has_bundle(bundle)
            if not ok:
                print(f"assemble.py: warning: no bundle found in {bundle}", file=sys.stderr)
        else:
            bundle = workdir / pid
            ok = fetch_bundle(project["repo"], bundle) and has_bundle(bundle)
        if not ok:
            _clean(site, pid)
            skipped.append(pid)
            continue
        data = assemble_project(site, project, bundle, skip, max_minor, drop)
        done.append(pid)
        print(f"assemble.py: {pid}: {data['repo']}@{data['sha'][:8] or '?'}")
    if skipped:
        print(f"assemble.py: skipped (no bundle): {', '.join(skipped)}", file=sys.stderr)
    if not done:
        # Not an error: before the first source repo publishes a bundle the site
        # is just the intro and quick start.
        print("::warning::assemble.py: no project bundles were assembled", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
