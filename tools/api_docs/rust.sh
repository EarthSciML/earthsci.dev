#!/usr/bin/env bash
# API reference, Rust half (tools/api_docs/README.md §1):
#   1. extract docs for every Rust spelling in the manifest from rustdoc's
#      JSON output -> build/api-docs/rust.json
#   2. build and run every example program in the Rust examples directory
#
# Usage: rust.sh [--repo-root DIR] [--extract-only | --examples-only]
#
# Reads [bindings.rust] from docs/api/config.toml: `path` (crate directory),
# `crate` (package name), and optionally `features` (enabled for rustdoc),
# `docs_rs` (docs.rs root URL), `docs_rs_features`, `nightly` (pinned rustdoc
# toolchain) and `examples` (examples package directory). Several crates are
# listed as [[bindings.rust.crates]] (each with those same per-crate keys);
# each is run through rustdoc and the results are merged into rust.json.
#
# Rustdoc's JSON output is a nightly-only flag. The toolchain that runs it is
# $RUSTDOC_TOOLCHAIN, default the pinned nightly (config `nightly`, default
# below); CI installs exactly that one. Locally, when the pinned nightly is not installed and
# RUSTDOC_TOOLCHAIN is unset, the script falls back to the default stable
# toolchain with RUSTC_BOOTSTRAP=1, which unlocks the same flag.
#
# The examples build with the default toolchain, or $EXAMPLES_TOOLCHAIN.
#
# All build output (cargo target dirs, temporary files) goes under build/.
set -euo pipefail

source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
parse_args "$@"
cd "$ROOT"

mode="$MODE"
PINNED_NIGHTLY="$(cfg bindings.rust.nightly --default nightly-2026-10-01)"
# One entry per crate; the first is the default crate (see README.md).
mapfile -t CRATE_DIRS < <("$PY" "$TOOLS/config.py" --repo-root "$ROOT" crates --field path)
mapfile -t CRATE_MODULES < <("$PY" "$TOOLS/config.py" --repo-root "$ROOT" crates --field module)
mapfile -t CRATE_DOCS_RS < <("$PY" "$TOOLS/config.py" --repo-root "$ROOT" crates --field docs_rs)
mapfile -t CRATE_DOCS_RS_FEATURES < <("$PY" "$TOOLS/config.py" --repo-root "$ROOT" crates --field docs_rs_features)
mapfile -t CRATE_FEATURES < <("$PY" "$TOOLS/config.py" --repo-root "$ROOT" crates --field features)
ALL_MODULES="$(IFS=,; echo "${CRATE_MODULES[*]}")"

BUILD="$ROOT/build"
mkdir -p "$BUILD/api-docs" "$BUILD/tmp"
export TMPDIR="$BUILD/tmp"
export CARGO_INCREMENTAL=0

have_toolchain() {
    command -v rustup >/dev/null 2>&1 && rustup run "$1" rustc --version >/dev/null 2>&1
}

# `cargo` with an optional `+toolchain`, as an array so an empty one vanishes.
cargo_cmd() {
    if [[ -n "$1" ]]; then
        echo "cargo +$1"
    else
        echo "cargo"
    fi
}

extract() {
    local tc="${RUSTDOC_TOOLCHAIN:-$PINNED_NIGHTLY}"
    local bootstrap=0
    if ! have_toolchain "$tc"; then
        if [[ -n "${RUSTDOC_TOOLCHAIN:-}" ]]; then
            echo "rust.sh: toolchain '$tc' is not installed (rustup toolchain install $tc --profile minimal)" >&2
            exit 1
        fi
        echo "rust.sh: $PINNED_NIGHTLY is not installed; using the default toolchain with RUSTC_BOOTSTRAP=1" >&2
        tc=""
        bootstrap=1
    elif [[ "$tc" != nightly* ]]; then
        bootstrap=1
    fi

    local cargo rustdoc_version
    cargo="$(cargo_cmd "$tc")"
    if [[ -n "$tc" ]]; then
        rustdoc_version="$(rustup run "$tc" rustdoc --version)"
    else
        rustdoc_version="$(rustdoc --version)"
    fi
    echo "rust.sh: extracting with $rustdoc_version" >&2

    local target="$BUILD/api-docs/rust-rustdoc-target"
    if [[ "$bootstrap" == 1 ]]; then export RUSTC_BOOTSTRAP=1; fi
    local note=""
    if [[ "$bootstrap" == 1 ]]; then note=", RUSTC_BOOTSTRAP=1"; fi

    local partials=() i
    for i in "${!CRATE_DIRS[@]}"; do
        local dir="${CRATE_DIRS[$i]}" module="${CRATE_MODULES[$i]}"
        # `features` gate root exports named in the manifest; building with them
        # documents every listed spelling.
        local features="${CRATE_FEATURES[$i]}"
        # Pin to Cargo.lock when the crate commits one; libraries often do not.
        local locked=""
        if [[ -f "$dir/Cargo.lock" ]]; then locked="--locked"; fi
        CARGO_TARGET_DIR="$target" $cargo rustdoc \
            --manifest-path "$dir/Cargo.toml" \
            --lib $locked ${features:+--features "$features"} \
            -- -Z unstable-options --output-format json
        local out="$BUILD/api-docs/rust.json"
        if (( ${#CRATE_DIRS[@]} > 1 )); then
            out="$BUILD/api-docs/rust-$module.partial.json"
            partials+=("$out")
        fi
        "$PY" "$TOOLS/rust_extract.py" \
            --rustdoc-json "$target/doc/$module.json" \
            --manifest "$(cfg manifest)" \
            --cargo-toml "$dir/Cargo.toml" \
            --crate-dir "$dir" --crate-module "$module" --crates "$ALL_MODULES" \
            --docs-rs-root "${CRATE_DOCS_RS[$i]}" --docs-rs-features "${CRATE_DOCS_RS_FEATURES[$i]}" \
            --extractor "rustdoc JSON ($rustdoc_version$note)" \
            --out "$out"
    done
    unset RUSTC_BOOTSTRAP
    if (( ${#partials[@]} )); then
        "$PY" "$TOOLS/rust_extract.py" --merge "${partials[@]}" --manifest "$(cfg manifest)" \
            --out "$BUILD/api-docs/rust.json"
        rm -f "${partials[@]}"
    fi
}

examples() {
    local dir
    dir="$(examples_dir rust)"
    local cargo
    cargo="$(cargo_cmd "${EXAMPLES_TOOLCHAIN:-}")"
    local target="$BUILD/api-docs/rust-examples-target"

    # Every program in the directory must be declared as an [[example]], or
    # it would silently never run.
    local declared names=() missing=0
    declared="$(sed -n 's/^path = "\(.*\.rs\)"$/\1/p' "$dir/Cargo.toml")"
    for f in "$dir"/*.rs; do
        local base
        base="$(basename "$f")"
        [[ "$base" == support.rs ]] && continue
        if ! grep -qx "$base" <<<"$declared"; then
            echo "rust.sh: $dir/$base is not declared as an [[example]] in $dir/Cargo.toml" >&2
            missing=1
        fi
        names+=("${base%.rs}")
    done
    [[ "$missing" == 0 ]] || exit 1

    # The crate's own lock file pins the examples' dependencies, so a
    # dependency change in the crate never leaves a second lock file stale.
    # Cargo then adds only the examples package itself (the copy is gitignored).
    # A crate inside a workspace has its lock file at the workspace root, so
    # look upward from the crate directory.
    local lock="${CRATE_DIRS[0]}" ; lock="${lock%/}"
    while [[ ! -f "$lock/Cargo.lock" && "$lock" == */* ]]; do lock="${lock%/*}"; done
    [[ -f "$lock/Cargo.lock" ]] || { echo "rust.sh: no Cargo.lock at or above ${CRATE_DIRS[0]}" >&2; exit 1; }
    cp "$lock/Cargo.lock" "$dir/Cargo.lock"
    CARGO_TARGET_DIR="$target" $cargo build --manifest-path "$dir/Cargo.toml" --examples

    # Run from an unrelated directory: the examples find their fixtures
    # relative to their own source, never the working directory.
    local failed=()
    for name in "${names[@]}"; do
        echo "--- example: $name" >&2
        if ! (cd "$TMPDIR" && "$target/debug/examples/$name"); then
            failed+=("$name")
        fi
    done
    if (( ${#failed[@]} )); then
        echo "rust.sh: ${#failed[@]} of ${#names[@]} examples failed: ${failed[*]}" >&2
        exit 1
    fi
    echo "rust.sh: all ${#names[@]} examples passed" >&2
}

if [[ "$mode" != examples ]]; then extract; fi
if [[ "$mode" != extract ]]; then examples; fi
