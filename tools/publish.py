#!/usr/bin/env python3
"""Merge a freshly built docs bundle into the source repo's `docs-bundle` branch.

The branch is a single orphan commit holding the whole bundle: the latest
narrative content from main plus one rendered API tree per version
(tools/versions.py). Publishing a tag adds or replaces that tag's tree and keeps
everything else; publishing `dev` replaces the content and the dev tree. The
branch is rewritten as one commit, so history does not grow, but it is pushed
with a lease so two concurrent publishes (a main push and a tag push) cannot
overwrite each other: the loser re-clones, re-merges and retries.

    python3 tools/publish.py --bundle build/bundle --version dev --repo OWNER/NAME
    python3 tools/publish.py --bundle b --version v1.2.0 --existing DIR   # no network; merge into DIR only

Auth: the token in $GH_TOKEN (or $GITHUB_TOKEN).
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import versions  # noqa: E402

BRANCH = "docs-bundle"


def git(*args, cwd=None, check=True):
    return subprocess.run(["git", *args], cwd=cwd, check=check, capture_output=True, text=True,
                          env={**os.environ, "GIT_TERMINAL_PROMPT": "0"})


def clone_existing(url: str, dest: Path) -> str:
    """Clone the bundle branch into dest, drop .git; return its sha ('' if absent)."""
    r = git("clone", "--quiet", "--depth", "1", "--branch", BRANCH, url, str(dest), check=False)
    if r.returncode != 0:
        dest.mkdir(parents=True, exist_ok=True)
        return ""
    sha = git("rev-parse", "HEAD", cwd=dest).stdout.strip()
    shutil.rmtree(dest / ".git")
    return sha


def publish(bundle: Path, version: str, repo: str, token: str, max_minor: int,
            message: str, attempts: int = 5) -> dict:
    url = f"https://x-access-token:{token}@github.com/{repo}.git"
    for attempt in range(1, attempts + 1):
        work = Path(tempfile.mkdtemp(prefix="docs-bundle-"))
        try:
            old = clone_existing(url, work / "tree")
            info = versions.merge(work / "tree", bundle, version, max_minor)
            tree = work / "tree"
            git("init", "-q", "-b", BRANCH, cwd=tree)
            git("config", "user.name", "github-actions[bot]", cwd=tree)
            git("config", "user.email", "41898282+github-actions[bot]@users.noreply.github.com", cwd=tree)
            git("add", "-A", cwd=tree)
            git("commit", "-q", "-m", message, cwd=tree)
            r = git("push", "-q", f"--force-with-lease={BRANCH}:{old}", url, f"{BRANCH}:{BRANCH}",
                    cwd=tree, check=False)
            if r.returncode == 0:
                return info
            print(f"publish.py: push rejected (attempt {attempt}/{attempts}), retrying: "
                  f"{r.stderr.strip().splitlines()[-1] if r.stderr.strip() else ''}", file=sys.stderr)
        finally:
            shutil.rmtree(work, ignore_errors=True)
    raise SystemExit("publish.py: could not push the docs-bundle branch")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bundle", type=Path, required=True, help="bundle built by tools/bundle.py")
    ap.add_argument("--version", required=True, help="`dev` or a release tag vX.Y.Z")
    ap.add_argument("--repo", help="OWNER/NAME of the source repo")
    ap.add_argument("--existing", type=Path, help="merge into this directory instead of pushing")
    ap.add_argument("--max-minor-versions", type=int, default=8,
                    help="keep the newest patch of this many minor series (0 = all)")
    ap.add_argument("--message", default="docs bundle")
    args = ap.parse_args(argv)
    if args.existing:
        info = versions.merge(args.existing, args.bundle, args.version, args.max_minor_versions)
    else:
        token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
        if not (token and args.repo):
            ap.error("--repo and $GH_TOKEN are required unless --existing is given")
        info = publish(args.bundle, args.version, args.repo, token,
                       args.max_minor_versions, args.message)
    print(f"publish.py: {args.version} published; versions now: {', '.join(info['versions'])} "
          f"(latest: {info['latest'] or 'none'})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
