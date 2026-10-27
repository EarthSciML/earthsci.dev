#!/usr/bin/env python3
"""Extract tagged regions from a repository's example programs.

The examples live under the directory named by docs/api/config.toml (default
docs/api/examples). The region syntax is defined in tools/api_docs/README.md:

    # region: load ops=load_path,load_string
    doc = earthsci_ast.load_path(FIXTURES / "decay.esm")
    setup_only()  # hide
    # endregion

The comment token is chosen by file extension. A region id must be unique
within one binding, regions may not nest, and every name in ``ops=`` must be a
canonical name from the manifest (api-surface.json) that the binding has a spelling for.

This is the only implementation of these rules: every binding script calls it.

Usage:

    python3 regions.py --repo-root .                  # check every configured binding
    python3 regions.py --repo-root . --binding go     # check one
    python3 regions.py --repo-root . --json           # dump the regions

Exits non-zero, listing every problem, if any file breaks the rules. The
binding scripts may call it to check their examples; render.py imports it.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import config as config_mod  # noqa: E402

BINDINGS = config_mod.SUPPORTED  # every language with an extractor; a repo configures a subset

# Source extensions per binding, and the line-comment token for each extension.
EXTENSIONS = {
    "julia": (".jl",),
    "typescript": (".ts", ".mts", ".js", ".mjs"),
    "python": (".py",),
    "rust": (".rs",),
    "go": (".go",),
}
COMMENT_TOKEN = {
    ".jl": "#",
    ".py": "#",
    ".ts": "//",
    ".mts": "//",
    ".js": "//",
    ".mjs": "//",
    ".rs": "//",
    ".go": "//",
}
# Directories that hold build products or dependencies, never examples.
SKIP_DIRS = {"node_modules", "target", "build", "dist", "__pycache__", "vendor"}

ID_RE = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9_.-]*$")


@dataclass
class Region:
    binding: str
    id: str
    ops: list[str]
    path: str  # relative to the repository root
    start_line: int  # 1-based line of the first shown code line
    end_line: int  # 1-based line of the last shown code line
    code: str  # dedented, hide lines removed
    language: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class ExtractResult:
    regions: list[Region] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)


def _markers(token: str):
    t = re.escape(token)
    start = re.compile(rf"^\s*{t}\s*region:\s*(?P<rest>.*?)\s*$")
    end = re.compile(rf"^\s*{t}\s*endregion\b.*$")
    hide = re.compile(rf"\s*{t}\s*hide\s*$")
    return start, end, hide


def dedent(lines: list[str]) -> list[str]:
    """Remove the common leading whitespace of the non-blank lines."""
    widths = [len(ln) - len(ln.lstrip()) for ln in lines if ln.strip()]
    cut = min(widths) if widths else 0
    out = [ln[cut:] if ln.strip() else "" for ln in lines]
    while out and not out[0].strip():
        out.pop(0)
    while out and not out[-1].strip():
        out.pop()
    return out


def extract_file(
    binding: str,
    path: Path,
    rel: str,
    canonical: dict[str, set[str]] | None,
    result: ExtractResult,
) -> None:
    token = COMMENT_TOKEN[path.suffix]
    start_re, end_re, hide_re = _markers(token)
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as exc:
        result.errors.append(f"{rel}: cannot read: {exc}")
        return

    open_region: dict | None = None
    for lineno, raw in enumerate(text.splitlines(), start=1):
        m = start_re.match(raw)
        if m:
            if open_region is not None:
                result.errors.append(
                    f"{rel}:{lineno}: region {m.group('rest')!r} opens inside region "
                    f"{open_region['id']!r} (line {open_region['line']}); regions may not nest"
                )
                continue
            rid, ops, problems = parse_header(m.group("rest"), canonical, binding)
            for p in problems:
                result.errors.append(f"{rel}:{lineno}: {p}")
            open_region = {"id": rid, "ops": ops, "line": lineno, "body": []}
            continue
        if end_re.match(raw):
            if open_region is None:
                result.errors.append(f"{rel}:{lineno}: endregion without a region")
                continue
            body = open_region["body"]
            shown = [(n, ln) for (n, ln) in body if not hide_re.search(ln)]
            code = "\n".join(dedent([ln for (_, ln) in shown]))
            # The link range covers what the page shows, not hidden or blank lines.
            visible = [n for (n, ln) in shown if ln.strip()]
            first = visible[0] if visible else lineno
            last = visible[-1] if visible else lineno
            if not code.strip():
                result.errors.append(
                    f"{rel}:{open_region['line']}: region {open_region['id']!r} shows no code"
                )
            if open_region["id"]:
                result.regions.append(
                    Region(
                        binding=binding,
                        id=open_region["id"],
                        ops=open_region["ops"],
                        path=rel,
                        start_line=first,
                        end_line=last,
                        code=code,
                        language=binding,
                    )
                )
            open_region = None
            continue
        if open_region is not None:
            open_region["body"].append((lineno, raw.rstrip()))
    if open_region is not None:
        result.errors.append(
            f"{rel}:{open_region['line']}: region {open_region['id']!r} is never closed"
        )


def parse_header(
    rest: str, canonical: dict[str, set[str]] | None, binding: str
) -> tuple[str, list[str], list[str]]:
    """Parse ``<id> [ops=a,b]`` into (id, ops, problems)."""
    problems: list[str] = []
    parts = rest.split()
    if not parts:
        return "", [], ["region marker has no id"]
    rid = parts[0]
    if rid.startswith("ops="):
        return "", [], ["region marker has no id before ops="]
    if not ID_RE.match(rid):
        problems.append(f"region id {rid!r} has characters outside [A-Za-z0-9_.-]")
    ops: list[str] = []
    for extra in parts[1:]:
        if extra.startswith("ops="):
            ops.extend(o for o in extra[4:].split(",") if o)
        else:
            problems.append(f"region {rid!r}: unexpected token {extra!r} (only ops= is allowed)")
    if canonical is not None:
        for op in ops:
            if op not in canonical:
                problems.append(
                    f"region {rid!r}: ops= names {op!r}, which is not a canonical name "
                    "in the manifest (use the canonical name, not the binding spelling)"
                )
            elif binding not in canonical[op]:
                problems.append(
                    f"region {rid!r}: ops= names {op!r}, which has no {binding} spelling"
                )
    return rid, ops, problems


def example_files(examples_dir: Path, binding: str, base: Path | None = None) -> list[Path]:
    base = base or examples_dir / binding
    if not base.is_dir():
        return []
    exts = EXTENSIONS[binding]
    files: list[Path] = []
    for p in sorted(base.rglob("*")):
        rel_parts = p.relative_to(base).parts
        if any(part in SKIP_DIRS or part.startswith(".") for part in rel_parts[:-1]):
            continue
        if p.is_file() and p.suffix in exts:
            files.append(p)
    return files


def canonical_names(root: Path, manifest: str = "api-surface.json") -> dict[str, set[str]]:
    """Map each canonical name in the manifest to the bindings that spell it."""
    surface = json.loads((root / manifest).read_text(encoding="utf-8"))
    names: dict[str, set[str]] = {}
    for s in surface["symbols"]:
        names.setdefault(s["name"], set()).update(s["bindings"])
    return names


def extract(
    root: Path,
    bindings=BINDINGS,
    canonical: dict[str, set[str]] | None = None,
    examples_dir: Path | None = None,
    manifest: str = "api-surface.json",
    binding_dirs: dict[str, Path] | None = None,
) -> ExtractResult:
    """Extract every region of the given bindings, collecting all problems.

    ``examples_dir`` defaults to ``<root>/docs/api/examples``; region paths are
    reported relative to ``root``. ``binding_dirs`` overrides the directory of
    individual bindings (config ``bindings.<name>.examples``).
    """
    if examples_dir is None:
        examples_dir = root / "docs" / "api" / "examples"
    if canonical is None:
        canonical = canonical_names(root, manifest)
    result = ExtractResult()
    for binding in bindings:
        before = len(result.regions)
        for path in example_files(examples_dir, binding, (binding_dirs or {}).get(binding)):
            try:
                rel = path.resolve().relative_to(root.resolve()).as_posix()
            except ValueError:
                rel = path.as_posix()
            extract_file(binding, path, rel, canonical, result)
        seen: dict[str, Region] = {}
        for region in result.regions[before:]:
            prior = seen.get(region.id)
            if prior is not None:
                result.errors.append(
                    f"{region.path}:{region.start_line}: region id {region.id!r} already used "
                    f"by {prior.path}:{prior.start_line} in the {binding} examples"
                )
            else:
                seen[region.id] = region
    return result


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "--repo-root", "--root", dest="root", type=Path, default=Path.cwd(),
        help="repository root (default: cwd)",
    )
    ap.add_argument(
        "--examples-dir", type=Path, help="examples directory (default: from docs/api/config.toml)"
    )
    ap.add_argument(
        "--binding", choices=BINDINGS, action="append", help="limit to one binding (repeatable)"
    )
    ap.add_argument("--json", action="store_true", help="print the extracted regions as JSON")
    args = ap.parse_args(argv)

    try:
        cfg = config_mod.load(args.root)
    except config_mod.ConfigError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    bindings = tuple(args.binding or cfg.names)
    result = extract(
        cfg.root,
        bindings,
        examples_dir=args.examples_dir or cfg.path(cfg.examples),
        manifest=cfg.manifest,
        binding_dirs={b: cfg.binding_dir(b) for b in cfg.names},
    )
    if args.json:
        json.dump([r.to_dict() for r in result.regions], sys.stdout, indent=2)
        sys.stdout.write("\n")
    else:
        for r in result.regions:
            print(f"{r.binding:10} {r.id:20} {r.path}:{r.start_line}  ops={','.join(r.ops)}")
    for err in result.errors:
        print(f"error: {err}", file=sys.stderr)
    return 1 if result.errors else 0


if __name__ == "__main__":
    sys.exit(main())
