#!/usr/bin/env python3
"""The per-repository API-doc configuration, ``<repo>/docs/api/config.toml``.

Every tool in this directory reads its project-specific settings here, so
nothing about a particular repository is hard-coded in the tooling. The schema
is documented in README.md. This module is also a small CLI for shell scripts
and CI:

    config.py --repo-root . get bindings.python.path
    config.py --repo-root . get bindings.rust.features --join ,
    config.py --repo-root . bindings            # binding names, one per line
    config.py --repo-root . bindings --json     # JSON array, for a job matrix

Python 3.11 or newer (tomllib), standard library only.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:  # Python < 3.11
    sys.exit("the API-doc tooling needs Python 3.11 or newer (for tomllib)")

# The languages the tooling has extractors for, with their display names, the
# source extensions of their example programs, and the comment token of each.
SUPPORTED = ("julia", "typescript", "python", "rust", "go")
LABELS = {
    "julia": "Julia",
    "typescript": "TypeScript",
    "python": "Python",
    "rust": "Rust",
    "go": "Go",
}


_CRATE_KEYS = ("path", "crate", "features", "docs_rs", "docs_rs_features")


class ConfigError(Exception):
    pass


@dataclass
class Config:
    root: Path
    manifest: str = "api-surface.json"
    examples: str = "docs/api/examples"
    groups: str = "docs/api/groups.toml"
    guides: str = "docs/api/guides"
    out: str = "docs/content/generated/api"
    reference: dict = field(default_factory=dict)
    bindings: dict[str, dict] = field(default_factory=dict)
    site: dict = field(default_factory=dict)

    # ---- derived values --------------------------------------------------

    @property
    def names(self) -> list[str]:
        """Configured binding names, in configuration order."""
        return list(self.bindings)

    def label(self, binding: str) -> str:
        return self.bindings[binding].get("label") or LABELS[binding]

    @property
    def title(self) -> str:
        """The project's display name (``reference.title``, else ``site.toml``'s)."""
        return self.reference.get("title") or self.site.get("title") or self.root.resolve().name

    @property
    def repo_url(self) -> str:
        url = self.reference.get("repo_url")
        if url:
            return url.rstrip("/")
        repo = self.site.get("repo")
        if repo:
            return f"https://github.com/{repo}"
        return ""

    def path(self, rel: str) -> Path:
        return self.root / rel

    def binding_dir(self, binding: str) -> Path:
        """Where this binding's example programs live."""
        return self.root / self.bindings[binding].get(
            "examples", f"{self.examples}/{binding}"
        )

    def rust_crates(self) -> list[dict]:
        """The Rust crates to document, defaults filled in, first = default crate.

        ``[[bindings.rust.crates]]`` lists them; the legacy single-crate keys
        (``path``, ``crate``, ``features``, ``docs_rs``, ``docs_rs_features``
        directly under ``[bindings.rust]``) are one crate.
        """
        b = self.bindings.get("rust")
        if b is None:
            return []
        entries = b.get("crates") or [{k: b[k] for k in _CRATE_KEYS if k in b}]
        crates = []
        for e in entries:
            name = e.get("crate") or Path(e["path"]).name
            crates.append(
                {
                    "path": str(e["path"]).strip("/"),
                    "crate": name,
                    "module": name.replace("-", "_"),
                    "features": list(e.get("features", [])),
                    "docs_rs": e.get("docs_rs") or f"https://docs.rs/{name}/latest/",
                    "docs_rs_features": list(e.get("docs_rs_features", [])),
                }
            )
        return crates

    def get(self, dotted: str):
        """A value by dotted key, e.g. ``bindings.rust.crate``; None when absent."""
        node = {
            "manifest": self.manifest,
            "examples": self.examples,
            "groups": self.groups,
            "guides": self.guides,
            "out": self.out,
            "reference": self.reference,
            "bindings": self.bindings,
        }
        for part in dotted.split("."):
            if not isinstance(node, dict) or part not in node:
                return None
            node = node[part]
        return node


def config_path(root: Path) -> Path:
    return root / "docs" / "api" / "config.toml"


def load(root: Path) -> Config:
    root = Path(root).resolve()
    path = config_path(root)
    if not path.is_file():
        raise ConfigError(f"{path}: not found (see tools/api_docs/README.md for the schema)")
    try:
        data = tomllib.loads(path.read_text(encoding="utf-8"))
    except tomllib.TOMLDecodeError as exc:
        raise ConfigError(f"{path}: {exc}") from exc
    cfg = Config(root=root)
    for key in ("manifest", "examples", "groups", "guides", "out"):
        if key in data:
            setattr(cfg, key, str(data[key]))
    cfg.reference = dict(data.get("reference", {}))
    cfg.bindings = {k: dict(v) for k, v in data.get("bindings", {}).items()}
    if not cfg.bindings:
        raise ConfigError(f"{path}: no [bindings.<name>] tables")
    for name, b in cfg.bindings.items():
        if name not in SUPPORTED:
            raise ConfigError(f"{path}: unknown binding {name!r} (supported: {', '.join(SUPPORTED)})")
        if name == "rust" and "crates" in b:
            if not isinstance(b["crates"], list) or not b["crates"]:
                raise ConfigError(f"{path}: [bindings.rust] `crates` must be a non-empty [[bindings.rust.crates]] list")
            if "path" in b:
                raise ConfigError(f"{path}: [bindings.rust] has both `path` and `crates`; use one form")
            for c in b["crates"]:
                if "path" not in c:
                    raise ConfigError(f"{path}: every [[bindings.rust.crates]] needs `path`")
            mods = [c.get("crate") or Path(c["path"]).name for c in b["crates"]]
            if len(set(m.replace("-", "_") for m in mods)) != len(mods):
                raise ConfigError(f"{path}: duplicate crate in [[bindings.rust.crates]]")
        elif "path" not in b:
            raise ConfigError(f"{path}: [bindings.{name}] needs `path`")
    site_toml = root / "docs" / "site.toml"
    if site_toml.is_file():
        try:
            cfg.site = tomllib.loads(site_toml.read_text(encoding="utf-8"))
        except tomllib.TOMLDecodeError as exc:
            raise ConfigError(f"{site_toml}: {exc}") from exc
    return cfg


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo-root", type=Path, default=Path.cwd())
    sub = ap.add_subparsers(dest="cmd", required=True)
    g = sub.add_parser("get", help="print one value")
    g.add_argument("key")
    g.add_argument("--default", default=None)
    g.add_argument("--join", default="\n", help="separator for list values")
    b = sub.add_parser("bindings", help="list configured binding names")
    b.add_argument("--json", action="store_true")
    c = sub.add_parser("crates", help="the Rust crates, one per line (or JSON)")
    c.add_argument("--field", default="path", help="crate, module, path, docs_rs, features, docs_rs_features")
    c.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)
    try:
        cfg = load(args.repo_root)
    except ConfigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    if args.cmd == "bindings":
        if args.json:
            print(json.dumps(cfg.names))
        else:
            print("\n".join(cfg.names))
        return 0
    if args.cmd == "crates":
        crates = cfg.rust_crates()
        if args.json:
            print(json.dumps(crates))
        else:
            for cr in crates:
                v = cr[args.field]
                print(",".join(v) if isinstance(v, list) else v)
        return 0
    val = cfg.get(args.key)
    if val is None:
        if args.default is None:
            print(f"error: no config value {args.key!r}", file=sys.stderr)
            return 1
        val = args.default
    if isinstance(val, list):
        val = args.join.join(str(v) for v in val)
    elif isinstance(val, bool):
        val = "true" if val else "false"
    print(val)
    return 0


if __name__ == "__main__":
    sys.exit(main())
