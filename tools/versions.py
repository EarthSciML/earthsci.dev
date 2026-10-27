#!/usr/bin/env python3
"""Version bookkeeping for the versioned API reference (DESIGN.md, "Versions").

A docs bundle keeps one rendered API tree per version under `api/<version>/`:
`dev` (the main branch) and one `vX.Y.Z` per released tag. The set of versions
is whatever directories exist, so there is no separate state to drift. This
module orders and prunes them, merges a freshly built bundle into an existing
one (what the source repo's `docs-bundle` branch holds), and writes
`versions.json`:

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


def scan(api_dir: Path) -> list[str]:
    if not api_dir.is_dir():
        return []
    return order(p.name for p in api_dir.iterdir() if p.is_dir() and is_version(p.name))


def effective_latest(names) -> str | None:
    """The version served at /<id>/generated/api/: the newest release, falling
    back to `dev` for a project with no release yet."""
    return latest(names) or (DEV if DEV in set(names) else None)


def write_versions(bundle: Path, max_minor: int = 0) -> dict:
    """Prune `bundle/api/` to the versions worth keeping and write versions.json."""
    api = bundle / "api"
    found = scan(api)
    keep = set(prune(found, max_minor))
    for name in found:
        if name not in keep:
            shutil.rmtree(api / name)
    names = order(keep)
    info = {"latest": latest(names), "versions": names}
    if names:
        (bundle / "versions.json").write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    elif (bundle / "versions.json").exists():
        (bundle / "versions.json").unlink()
    return info


def _rm(p: Path) -> None:
    if p.is_dir():
        shutil.rmtree(p)
    elif p.exists():
        p.unlink()


def merge(existing: Path, new: Path, version: str, max_minor: int = 0) -> dict:
    """Merge the bundle just built for `version` into `existing` (in place).

    A `dev` build carries the whole bundle (bundle.json, content/, static/ and
    api/dev/): it replaces those. A tag build carries only api/<tag>/, so every
    other version, and the narrative content from main, stays as it was.
    """
    if not is_version(version):
        raise ValueError(f"{version!r} is not `dev` or a vX.Y.Z tag")
    existing.mkdir(parents=True, exist_ok=True)
    if version == DEV:
        for sub in ("content", "static", "bundle.json"):
            _rm(existing / sub)
            if (new / sub).exists():
                (shutil.copytree if (new / sub).is_dir() else shutil.copy2)(new / sub, existing / sub)
    elif not (existing / "bundle.json").exists() and (new / "bundle.json").exists():
        shutil.copy2(new / "bundle.json", existing / "bundle.json")
    _rm(existing / "api" / version)
    if (new / "api" / version).is_dir():
        shutil.copytree(new / "api" / version, existing / "api" / version, ignore=IGNORED)
    return write_versions(existing, max_minor)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bundle", type=Path, required=True)
    args = ap.parse_args(argv)
    names = scan(args.bundle / "api")
    print(json.dumps({"latest": latest(names), "versions": names}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
