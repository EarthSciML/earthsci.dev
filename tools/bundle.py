#!/usr/bin/env python3
"""Build a docs bundle from a source repo (DESIGN.md, "Source-repo docs contract").

A bundle is the complete docs of ONE version, a directory

    bundle.json   id, title, description, repo, nav, sha, version
    content/      the repo's Hugo content tree, including generated/ pages and
                  the API reference (generated/api)
    static/       the repo's static assets (optional)

`--version dev` (the default) is main; `--version vX.Y.Z` is a release tag, built
from the code and docs config at that tag. tools/publish.py merges the bundle
into the `docs-bundle` branch of the source repo as the directory `<version>/`
(tools/versions.py), and tools/assemble.py later copies every version into the
site.

    python3 tools/bundle.py --repo-root . --out build/bundle
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import versions  # noqa: E402

try:
    import tomllib
except ImportError:  # Python < 3.11
    import tomli as tomllib  # type: ignore

REQUIRED = ("id", "title", "description", "repo")
ID_RE = re.compile(r"^[a-z][a-z0-9-]*$")
IGNORED = shutil.ignore_patterns(".git", "__pycache__", ".DS_Store", "*.pyc")


class BundleError(Exception):
    pass


def load_site_toml(repo_root: Path) -> dict:
    path = repo_root / "docs" / "site.toml"
    if not path.is_file():
        raise BundleError(f"{path} not found (every docs repo needs docs/site.toml)")
    with path.open("rb") as f:
        cfg = tomllib.load(f)
    missing = [k for k in REQUIRED if not cfg.get(k)]
    if missing:
        raise BundleError(f"{path}: missing required keys: {', '.join(missing)}")
    if not ID_RE.match(cfg["id"]):
        raise BundleError(f"{path}: id {cfg['id']!r} must match {ID_RE.pattern}")
    for i, item in enumerate(cfg.get("nav", [])):
        if not isinstance(item, dict) or not item.get("title") or not item.get("path"):
            raise BundleError(f"{path}: nav[{i}] needs `title` and `path`")
    return cfg


def git_sha(repo_root: Path) -> str:
    try:
        return subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return os.environ.get("GITHUB_SHA", "")


def build_bundle(repo_root: Path, out: Path, version: str = versions.DEV) -> dict:
    if not versions.is_version(version):
        raise BundleError(f"--version {version!r} must be `dev` or a vX.Y.Z tag")
    repo_root = repo_root.resolve()
    cfg = load_site_toml(repo_root)
    content = repo_root / cfg.get("content", "docs/content")
    static = repo_root / cfg.get("static", "docs/static")
    if not content.is_dir():
        raise BundleError(f"content directory {content} not found")
    if not (content / "_index.md").is_file():
        print(f"bundle.py: warning: {content}/_index.md missing; the site will "
              "generate a bare landing page", file=sys.stderr)

    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    has_api = (content / "generated" / "api" / "_index.md").is_file()
    shutil.copytree(content, out / "content", ignore=IGNORED)
    if static.is_dir():
        shutil.copytree(static, out / "static", ignore=IGNORED)

    nav = [dict(n) for n in cfg.get("nav", [])]
    if has_api and not any(n["path"].strip("/") == "generated/api" for n in nav):
        nav.append({"title": "API", "path": "generated/api/"})
    elif not has_api and (repo_root / "docs" / "api" / "config.toml").is_file():
        print("bundle.py: warning: docs/api/config.toml exists but no generated "
              "API reference was found; run tools/api_docs/render.py first",
              file=sys.stderr)

    meta = {
        "id": cfg["id"],
        "title": cfg["title"],
        "description": cfg["description"],
        "repo": cfg["repo"],
        "nav": nav,
        "sha": git_sha(repo_root),
        "version": version,
    }
    (out / "bundle.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    return meta


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--repo-root", type=Path, default=Path("."))
    ap.add_argument("--out", type=Path, default=Path("build/bundle"))
    ap.add_argument("--version", default=versions.DEV,
                    help="`dev` (the default) or a release tag vX.Y.Z")
    args = ap.parse_args(argv)
    try:
        meta = build_bundle(args.repo_root, args.out, args.version)
    except BundleError as e:
        print(f"bundle.py: {e}", file=sys.stderr)
        return 1
    print(f"bundle.py: bundled {meta['id']} {args.version} ({meta['repo']}) -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
