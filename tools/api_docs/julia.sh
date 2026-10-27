#!/usr/bin/env bash
# Julia half of the generated API reference (tools/api_docs/README.md §1).
#
#   julia.sh [--repo-root DIR]                  extract docstrings, then run the examples
#   julia.sh --extract-only                     write build/api-docs/julia.json only
#   julia.sh --examples-only                    run the Julia examples only
#
# Both halves use the Julia environment named by [bindings.julia] `env` in
# docs/api/config.toml (default <examples>/julia). It must hold the package
# under documentation (typically dev-linked through `[sources]`, so Julia >= 1.11
# is required), JSON3, and any packages listed in `load` (extensions whose
# methods and docstrings are part of the surface). Set JULIA to choose the
# executable; JULIA_THREADS (default 2) sets --threads.
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
parse_args "$@"

ENV_DIR="$ROOT/$(cfg bindings.julia.env --default "$(cfg examples)/julia")"
PKG="$(cfg bindings.julia.package)"
LOAD="$(cfg bindings.julia.load --default '' --join ,)"
OUT_DIR="$ROOT/build/api-docs"
JULIA="${JULIA:-julia}"
JULIA_FLAGS=(--startup-file=no --threads="${JULIA_THREADS:-2}" --project="$ENV_DIR")

mkdir -p "$OUT_DIR"

echo "==> instantiating $ENV_DIR"
"$JULIA" "${JULIA_FLAGS[@]}" -e 'using Pkg; Pkg.instantiate()'

if [[ $MODE != examples ]]; then
  echo "==> extracting Julia docstrings"
  "$JULIA" "${JULIA_FLAGS[@]}" "$TOOLS/extract_julia.jl" \
    "$ROOT" "$PKG" "$ROOT/$(cfg manifest)" "$OUT_DIR/julia.json" "$LOAD"
fi

if [[ $MODE != extract ]]; then
  echo "==> running Julia examples"
  (cd "$ROOT" && "$JULIA" "${JULIA_FLAGS[@]}" "$TOOLS/run_julia_examples.jl" "$(examples_dir julia)")
fi
