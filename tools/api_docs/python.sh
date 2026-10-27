#!/usr/bin/env bash
# Python half of the generated API reference (tools/api_docs/README.md §1).
#
#   python.sh [--repo-root DIR]    extract docs, then run examples
#   python.sh --extract-only       write build/api-docs/python.json
#   python.sh --examples-only      run the Python examples (<examples>/python/*.py)
#
# Reads [bindings.python] from docs/api/config.toml: `path` (package
# directory, holding pyproject.toml), `package` (import name) and optionally
# `extras` (extras the examples need). A virtual environment under
# build/api-docs/venv-python holds the pinned extractor and the in-repo package
# (editable, with those extras). It is rebuilt when the package's
# pyproject.toml, the extras or the pinned griffe version change.
#
# Environment:
#   PYTHON           interpreter used to create the venv (default: the first of
#                    python3.13, python3.12, python3.11, python3 that is at
#                    least 3.11, which the pinned griffe requires)
#   API_DOCS_VENV    venv location (default: build/api-docs/venv-python)
set -euo pipefail

GRIFFE_VERSION="2.3.2"

source "$(dirname "${BASH_SOURCE[0]}")/lib.sh"
parse_args "$@"
cd "$ROOT"

PKG_DIR="$(cfg bindings.python.path)"
EXTRAS="$(cfg bindings.python.extras --default '' --join ,)"
mode="$MODE"

mkdir -p build/api-docs build/tmp
export TMPDIR="${TMPDIR:-$ROOT/build/tmp}"

VENV="${API_DOCS_VENV:-$ROOT/build/api-docs/venv-python}"
VPY="$VENV/bin/python"

stamp="griffe==$GRIFFE_VERSION [$EXTRAS] $(sha256sum "$PKG_DIR/pyproject.toml" | cut -d' ' -f1)"
if [[ ! -x "$VPY" || "$(cat "$VENV/.api-docs-stamp" 2>/dev/null)" != "$stamp" ]]; then
  echo "python.sh: creating $VENV" >&2
  "$PYTHON" -m venv --clear "$VENV"
  "$VPY" -m pip install --quiet --upgrade pip
  "$VPY" -m pip install --quiet "griffe==$GRIFFE_VERSION" -e "${PKG_DIR}${EXTRAS:+[$EXTRAS]}"
  echo "$stamp" > "$VENV/.api-docs-stamp"
fi

if [[ "$mode" != "examples" ]]; then
  "$VPY" "$TOOLS/extract_python.py" --repo-root "$ROOT" --output "$ROOT/build/api-docs/python.json"
fi

if [[ "$mode" != "extract" ]]; then
  # Run from a scratch directory so an example that found its fixtures through
  # the working directory, rather than its own location, fails here.
  rundir="$ROOT/build/api-docs/python-examples-run"
  rm -rf "$rundir"
  mkdir -p "$rundir"
  failed=()
  shopt -s nullglob
  EX_DIR="$(examples_dir python)"
  examples=("$EX_DIR"/*.py)
  if [[ ${#examples[@]} -eq 0 ]]; then
    echo "python.sh: no examples under $EX_DIR" >&2
    exit 1
  fi
  for ex in "${examples[@]}"; do
    name="$(basename "$ex")"
    if (cd "$rundir" && "$VPY" -W "error::DeprecationWarning:__main__" "$ex") > "$rundir/${name%.py}.log" 2>&1; then
      echo "ok    $name"
    else
      echo "FAIL  $name"
      sed 's/^/      /' "$rundir/${name%.py}.log"
      failed+=("$name")
    fi
  done
  if [[ ${#failed[@]} -gt 0 ]]; then
    echo "python.sh: ${#failed[@]} of ${#examples[@]} examples failed: ${failed[*]}" >&2
    exit 1
  fi
  echo "python.sh: all ${#examples[@]} examples passed"
fi
