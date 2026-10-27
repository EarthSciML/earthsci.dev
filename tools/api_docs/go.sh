#!/usr/bin/env bash
# Go half of the generated API reference (tools/api_docs/README.md §1).
#
#   go.sh [--repo-root DIR]     extract docstrings, then run examples
#   go.sh --extract-only        write build/api-docs/go.json only
#   go.sh --examples-only       run the Go examples only
#
# Reads [bindings.go] from docs/api/config.toml: `path` (the package directory
# to document) and `package` (its import path). Needs a Go toolchain and
# Python >= 3.11 for the region check. Works from any directory; all output
# goes under <repo>/build/.
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
parse_args "$@"

mkdir -p "$ROOT/build/api-docs" "$ROOT/build/tmp"
export TMPDIR="$ROOT/build/tmp"
export GOTMPDIR="$TMPDIR"

if [[ $MODE != examples ]]; then
  go -C "$TOOLS/extract_go" run . -root "$ROOT" -out "$ROOT/build/api-docs/go.json" \
    -manifest "$(cfg manifest)" -pkgdir "$(cfg bindings.go.path)" -import-path "$(cfg bindings.go.package)"
  echo "wrote build/api-docs/go.json"
fi

if [[ $MODE != extract ]]; then
  "$PY" "$TOOLS/regions.py" --repo-root "$ROOT" --binding go >/dev/null
  ex="$(examples_dir go)"
  go -C "$ex" vet ./...
  go -C "$ex" test -p 2 -count=1 ./...
fi
