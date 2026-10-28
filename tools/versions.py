#!/usr/bin/env python3
"""Version bookkeeping for versioned docs (DESIGN.md, "Versions").

A project's docs are published as one complete tree per version, in a directory
named for it: `dev` (the main branch) and one `vX.Y.Z` per released tag, each
holding that version's `bundle.json`, `content/` and `static/`. This is the
layout of the source repo's `docs-bundle` branch. The set of versions is
whatever such directories exist, so there is no separate state to drift. This
module orders and prunes them, merges a freshly built tree into an existing
bundle, migrates the older layout (one flat dev tree plus `api/<version>/`),
and writes `versions.json`:

    {"latest": "v1.2.0", "versions": ["dev", "v1.2.0", "v1.1.3"]}

`latest` is the newest non-prerelease tag, or null when there is none yet.

    python3 tools/versions.py --bundle DIR   # print the versions found in DIR
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import sys
from pathlib import Path

DEV = "dev"
_VERSION = re.compile(r"^v(\d+)\.(\d+)\.(\d+)(?:-([0-9A-Za-z.-]+))?$")
IGNORED = shutil.ignore_patterns(".git", "__pycache__", ".DS_Store", "*.pyc")


def parse(name: str):
    """(major, minor, patch, prerelease-or-None) for `vX.Y.Z[-pre]`, else None."""
    m = _VERSION.match(name)
    if not m:
        return None
    major, minor, patch, pre = m.groups()
    return int(major), int(minor), int(patch), pre


def is_version(name: str) -> bool:
    return name == DEV or parse(name) is not None


def is_release(name: str) -> bool:
    p = parse(name)
    return p is not None and p[3] is None


def _pre_key(pre):
    # semver: a release outranks its prereleases; numeric identifiers sort
    # numerically and below alphanumeric ones.
    if pre is None:
        return (1,)
    ids = [(0, int(i), "") if i.isdigit() else (1, 0, i) for i in pre.split(".")]
    return (0, tuple(ids))


def sort_key(name: str):
    p = parse(name)
    return (p[0], p[1], p[2], _pre_key(p[3]))


def order(names) -> list[str]:
    """`dev` first, then versions newest first; unrecognised names are dropped."""
    names = set(names)
    vers = sorted((n for n in names if parse(n)), key=sort_key, reverse=True)
    return ([DEV] if DEV in names else []) + vers


def latest(names) -> str | None:
    rel = sorted((n for n in names if is_release(n)), key=sort_key)
    return rel[-1] if rel else None


def prune(names, max_minor: int = 0) -> list[str]:
    """The versions worth keeping: the newest patch of each of the newest
    `max_minor` minor series (0 = no limit), plus any prerelease newer than the
    latest release. `dev` is always kept."""
    names = order(names)
    top = latest(names)
    rel = [n for n in names if is_release(n)]
    series: dict[tuple, str] = {}
    for n in rel:  # newest first, so the first per (major, minor) is its newest patch
        series.setdefault(parse(n)[:2], n)
    keep = list(series.values())
    if max_minor > 0:
        keep = keep[:max_minor]
    pres = [n for n in names if parse(n) and parse(n)[3] is not None
            and (top is None or sort_key(n) > sort_key(top))]
    return order([DEV] * (DEV in names) + keep + pres[:1])


def scan(root: Path) -> list[str]:
    """The versions published under `root`: directories named for a version that
    hold a bundle.json."""
    if not root.is_dir():
        return []
    return order(p.name for p in root.iterdir()
                 if p.is_dir() and is_version(p.name) and (p / "bundle.json").is_file())


def effective_latest(names) -> str | None:
    """The version served at /<id>/: the newest release, falling back to `dev`
    for a project with no release yet."""
    return latest(names) or (DEV if DEV in set(names) else None)


def write_versions(root: Path, max_minor: int = 0) -> dict:
    """Prune `root` to the versions worth keeping and write versions.json."""
    found = scan(root)
    keep = set(prune(found, max_minor))
    for name in found:
        if name not in keep:
            shutil.rmtree(root / name)
    names = order(keep)
    info = {"latest": latest(names), "versions": names}
    if names:
        (root / "versions.json").write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    elif (root / "versions.json").exists():
        (root / "versions.json").unlink()
    return info


def _rm(p: Path) -> None:
    if p.is_dir():
        shutil.rmtree(p)
    elif p.exists():
        p.unlink()


def is_legacy(root: Path) -> bool:
    """The first layout: one flat tree (bundle.json, content/, static/) for main
    plus api/<version>/ trees for the API reference alone."""
    return (root / "bundle.json").is_file() or (root / "api").is_dir()


def migrate_legacy(root: Path) -> None:
    """Rewrite a legacy bundle in place as a `dev` version tree. The flat tree
    was main's, and its api/dev tree becomes content/generated/api. The API-only
    trees of tags are dropped: they carry no narrative docs, so they are not
    versions of the whole project; the next tag build recreates them."""
    if (root / "bundle.json").is_file():
        tmp = root / ".migrate"
        _rm(tmp)
        tmp.mkdir()
        for sub in ("bundle.json", "content", "static"):
            if (root / sub).exists():
                shutil.move(str(root / sub), str(tmp / sub))
        if (root / "api" / "dev").is_dir():
            (tmp / "content" / "generated").mkdir(parents=True, exist_ok=True)
            _rm(tmp / "content" / "generated" / "api")
            shutil.copytree(root / "api" / "dev", tmp / "content" / "generated" / "api")
        _rm(root / DEV)
        tmp.rename(root / DEV)
    _rm(root / "api")
    _rm(root / "versions.json")


def merge(existing: Path, new: Path, version: str, max_minor: int = 0) -> dict:
    """Merge the tree just built for `version` into the bundle `existing` (in
    place): that version's directory is replaced, every other version stays."""
    if not is_version(version):
        raise ValueError(f"{version!r} is not `dev` or a vX.Y.Z tag")
    if not (new / "bundle.json").is_file():
        raise ValueError(f"{new} has no bundle.json (build it with tools/bundle.py)")
    existing.mkdir(parents=True, exist_ok=True)
    if is_legacy(existing):
        migrate_legacy(existing)
    _rm(existing / version)
    shutil.copytree(new, existing / version, ignore=IGNORED)
    return write_versions(existing, max_minor)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bundle", type=Path, required=True)
    args = ap.parse_args(argv)
    names = scan(args.bundle)
    print(json.dumps({"latest": latest(names), "versions": names}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
