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
# toolchain) and `examples` (examples package directory).
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
CRATE_DIR="$(cfg bindings.rust.path)"
CRATE="$(cfg bindings.rust.crate)"
CRATE_MODULE="${CRATE//-/_}"
FEATURES="$(cfg bindings.rust.features --default '' --join ,)"
DOCS_RS="$(cfg bindings.rust.docs_rs --default "https://docs.rs/$CRATE/latest/")"
DOCS_RS_FEATURES="$(cfg bindings.rust.docs_rs_features --default '' --join ,)"
PINNED_NIGHTLY="$(cfg bindings.rust.nightly --default nightly-2026-10-01)"

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
    # `features` gate root exports named in the manifest; building with them
    # documents every listed spelling.
    if [[ "$bootstrap" == 1 ]]; then export RUSTC_BOOTSTRAP=1; fi
    # Pin to Cargo.lock when the crate commits one; libraries often do not.
    local locked=""
    if [[ -f "$CRATE_DIR/Cargo.lock" ]]; then locked="--locked"; fi
    CARGO_TARGET_DIR="$target" $cargo rustdoc \
        --manifest-path "$CRATE_DIR/Cargo.toml" \
        --lib $locked ${FEATURES:+--features "$FEATURES"} \
        -- -Z unstable-options --output-format json
    unset RUSTC_BOOTSTRAP

    local note=""
    if [[ "$bootstrap" == 1 ]]; then note=", RUSTC_BOOTSTRAP=1"; fi
    "$PY" "$TOOLS/rust_extract.py" \
        --rustdoc-json "$target/doc/$CRATE_MODULE.json" \
        --manifest "$(cfg manifest)" \
        --cargo-toml "$CRATE_DIR/Cargo.toml" \
        --crate-dir "$CRATE_DIR" --crate-module "$CRATE_MODULE" \
        --docs-rs-root "$DOCS_RS" --docs-rs-features "$DOCS_RS_FEATURES" \
        --extractor "rustdoc JSON ($rustdoc_version$note)" \
        --out "$BUILD/api-docs/rust.json"
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
    local lock="$CRATE_DIR" ; lock="${lock%/}"
    while [[ ! -f "$lock/Cargo.lock" && "$lock" == */* ]]; do lock="${lock%/*}"; done
    [[ -f "$lock/Cargo.lock" ]] || { echo "rust.sh: no Cargo.lock at or above $CRATE_DIR" >&2; exit 1; }
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
