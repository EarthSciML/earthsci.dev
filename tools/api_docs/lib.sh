# Shared by <binding>.sh. Source it, then call `parse_args "$@"`.
#
# Sets TOOLS (this directory), ROOT (the repository being documented),
# MODE (all | extract | examples) and PY (a Python >= 3.11), and defines
# `cfg KEY [--join SEP]`, which reads a value from ROOT/docs/api/config.toml.

TOOLS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

parse_args() {
  ROOT="${REPO_ROOT:-$PWD}"
  MODE=all
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --repo-root) ROOT="$2"; shift 2 ;;
      --extract-only) MODE=extract; shift ;;
      --examples-only) MODE=examples; shift ;;
      -h|--help) echo "usage: $0 [--repo-root DIR] [--extract-only | --examples-only]"; exit 0 ;;
      *) echo "usage: $0 [--repo-root DIR] [--extract-only | --examples-only]" >&2; exit 2 ;;
    esac
  done
  ROOT="$(cd "$ROOT" && pwd)"
  find_python
}

# config.py, the region checker and the renderer need Python >= 3.11 (tomllib).
find_python() {
  if [[ -z "${PYTHON:-}" ]]; then
    # Plain python3 first: it is the one setup-python (or an active venv) put on PATH.
    local candidate
    for candidate in python3 python3.13 python3.12 python3.11; do
      if command -v "$candidate" >/dev/null 2>&1 &&
         "$candidate" -c 'import sys; sys.exit(sys.version_info < (3, 11))' 2>/dev/null; then
        PYTHON="$candidate"
        break
      fi
    done
  fi
  if [[ -z "${PYTHON:-}" ]]; then
    echo "$(basename "$0"): needs Python >= 3.11 (set PYTHON=/path/to/python3.11)" >&2
    exit 1
  fi
  PY="$PYTHON"
}

cfg() {
  "$PY" "$TOOLS/config.py" --repo-root "$ROOT" get "$@"
}

# The example directory of a binding, as an absolute path.
examples_dir() { # <binding>
  local custom
  custom="$(cfg "bindings.$1.examples" --default '')"
  if [[ -n "$custom" ]]; then echo "$ROOT/$custom"; else echo "$ROOT/$(cfg examples)/$1"; fi
}
