#!/usr/bin/env bash
# API reference, TypeScript binding (tools/api_docs/README.md §1):
#   1. extract docs for every TypeScript spelling in the manifest
#      -> build/api-docs/typescript.json
#   2. check the example region markers (regions.py), then type-check and run
#      every example under the TypeScript examples directory
#
# Usage: typescript.sh [--repo-root DIR] [--extract-only | --examples-only]
# Reads [bindings.typescript] from docs/api/config.toml: `path` (package
# directory), `package` (npm name), and optionally `entry` (default
# src/index.ts) and `dist` (built entry the examples import, default
# dist/esm/index.js). Needs Node >= 22.6 (examples run with Node's TypeScript
# type stripping).
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
parse_args "$@"

PKG_REL="$(cfg bindings.typescript.path)"
PKG="$ROOT/$PKG_REL"
NPM_NAME="$(cfg bindings.typescript.package)"
ENTRY="$(cfg bindings.typescript.entry --default src/index.ts)"
DIST="$(cfg bindings.typescript.dist --default dist/esm/index.js)"

mkdir -p "$ROOT/build/api-docs" "$ROOT/build/tmp"
export TMPDIR="$ROOT/build/tmp"

node -e 'const [a,b]=process.versions.node.split(".").map(Number); process.exit(a>22||(a===22&&b>=6)?0:1)' || {
  echo "typescript.sh: Node >= 22.6 required (found $(node --version))" >&2; exit 1; }

if [ ! -d "$PKG/node_modules" ]; then
  (cd "$PKG" && npm ci --no-audit --no-fund)
fi

if [ "$MODE" != "examples" ]; then
  node "$TOOLS/typescript_extract.mjs" "$ROOT" "$ROOT/build/api-docs/typescript.json" \
    "$PKG_REL" "$(cfg manifest)" "$ENTRY"
fi

if [ "$MODE" != "extract" ]; then
  # The examples import the built package; rebuild when the build is missing
  # or older than any source file.
  dist="$PKG/$DIST"
  if [ ! -f "$dist" ] || [ -n "$(find "$PKG/src" "$PKG/package.json" -newer "$dist" -print -quit)" ]; then
    (cd "$PKG" && npm run --silent build >"$ROOT/build/api-docs/typescript-build.log" 2>&1) || {
      cat "$ROOT/build/api-docs/typescript-build.log" >&2; exit 1; }
  fi
  "$PY" "$TOOLS/regions.py" --repo-root "$ROOT" --binding typescript >/dev/null
  node "$TOOLS/typescript_examples.mjs" "$ROOT" "$PKG_REL" "$(examples_dir typescript | sed "s#^$ROOT/##")" "$NPM_NAME" "$DIST"
fi
