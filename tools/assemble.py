#!/usr/bin/env python3
"""Assemble every project's docs bundle into this Hugo site.

For each [[project]] in projects.toml: obtain its bundle (the `docs-bundle`
branch of its repo, or a local directory given with --local), copy
content/ -> content/<id>/ and static/ -> static/<id>/, prefix root-relative
links with /<id>, and write data/projects/<id>.json for the header. A project
with no bundle is skipped with a warning.

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
from pathlib import Path

try:
    import tomllib
except ImportError:  # Python < 3.11
    import tomli as tomllib  # type: ignore

BRANCH = "docs-bundle"
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


def _prefixer(project_id: str, skip: set[str]):
    def repl(m: re.Match, text: str) -> str:
        rest = text[m.end():]
        first = re.split(r"[/?#)\"'\s]", rest, maxsplit=1)[0]
        if first == project_id or first in skip:
            return m.group(0)
        return f"{m.group(1)}/{project_id}/"
    return repl


def rewrite_links(text: str, project_id: str, skip: set[str] = frozenset()) -> str:
    """Prefix root-relative links with /<project_id>, leaving code alone."""
    repl = _prefixer(project_id, set(skip))
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
            for pat in _LINK_PATTERNS:
                seg = pat.sub(lambda mm, seg=seg: repl(mm, seg), seg)
            parts[i] = seg
        out.append("".join(parts))
    return "".join(out)


def rewrite_tree(root: Path, project_id: str, skip: set[str]) -> int:
    n = 0
    for p in root.rglob("*"):
        if p.suffix in TEXT_SUFFIXES and p.is_file():
            old = p.read_text(encoding="utf-8")
            new = rewrite_links(old, project_id, skip)
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


def assemble_project(site: Path, project: dict, bundle: Path, skip: set[str]) -> dict:
    pid = project["id"]
    meta = json.loads((bundle / "bundle.json").read_text(encoding="utf-8"))
    if meta["id"] != pid:
        raise SystemExit(f"assemble.py: bundle for {project['repo']} has id "
                         f"{meta['id']!r}, projects.toml says {pid!r}")
    _clean(site, pid)
    shutil.copytree(bundle / "content", site / "content" / pid)
    if (bundle / "static").is_dir():
        shutil.copytree(bundle / "static", site / "static" / pid)
    rewrite_tree(site / "content" / pid, pid, skip)

    index = site / "content" / pid / "_index.md"
    if not index.exists():
        index.write_text(
            f'---\ntitle: "{meta["title"]}"\ndescription: "{meta["description"]}"\n---\n',
            encoding="utf-8",
        )
    data = {
        "id": pid,
        "title": meta.get("title") or project["title"],
        "description": meta["description"],
        "repo": meta.get("repo") or project["repo"],
        "order": project.get("order", 100),
        "nav": meta.get("nav", []),
        "sha": meta.get("sha", ""),
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
    local = dict(item.split("=", 1) for item in args.local)
    workdir = args.workdir or site / "build" / "bundles"

    done, skipped = [], []
    for project in cfg["project"]:
        pid = project["id"]
        if args.only and pid not in args.only:
            continue
        if pid in local:
            bundle = Path(local[pid])
            ok = (bundle / "bundle.json").is_file()
            if not ok:
                print(f"assemble.py: warning: {bundle}/bundle.json not found", file=sys.stderr)
        else:
            bundle = workdir / pid
            ok = fetch_bundle(project["repo"], bundle) and (bundle / "bundle.json").is_file()
        if not ok:
            _clean(site, pid)
            skipped.append(pid)
            continue
        data = assemble_project(site, project, bundle, skip)
        done.append(pid)
        print(f"assemble.py: {pid}: {data['repo']}@{data['sha'][:8] or '?'}")
    if skipped:
        print(f"assemble.py: skipped (no bundle): {', '.join(skipped)}", file=sys.stderr)
    if not done:
        print("assemble.py: no project bundles were assembled", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
