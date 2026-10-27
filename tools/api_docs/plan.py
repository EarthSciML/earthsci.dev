#!/usr/bin/env python3
"""Print the configured bindings as a JSON array, for a job matrix.

    echo "bindings=$(python3 plan.py --repo-root .)" >> "$GITHUB_OUTPUT"

Prints ``[]`` when the repository has no ``docs/api/config.toml`` (no API
reference), so a caller's matrix is simply empty.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as config_mod  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo-root", type=Path, default=Path.cwd())
    args = ap.parse_args(argv)
    if not config_mod.config_path(args.repo_root.resolve()).is_file():
        print("[]")
        return 0
    try:
        cfg = config_mod.load(args.repo_root)
    except config_mod.ConfigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(cfg.names))
    return 0


if __name__ == "__main__":
    sys.exit(main())
