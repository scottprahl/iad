#!/bin/sh

set -eu

fail() {
    echo "FAIL: $*" >&2
    exit 1
}

case "${0}" in
    */*) CLI_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd) ;;
    *) CLI_DIR=$(pwd) ;;
esac

ROOT_DIR=${ROOT_DIR:-$(CDPATH= cd -- "$CLI_DIR/../.." && pwd)}
IAD_EXECUTABLE=${IAD_EXECUTABLE:-"$ROOT_DIR/iad"}
AD_EXECUTABLE=${AD_EXECUTABLE:-"$ROOT_DIR/ad"}
PYTHON=${PYTHON:-python3}

case "$IAD_EXECUTABLE" in
    /*) ;;
    *) IAD_EXECUTABLE="$ROOT_DIR/$IAD_EXECUTABLE" ;;
esac

case "$AD_EXECUTABLE" in
    /*) ;;
    *) AD_EXECUTABLE="$ROOT_DIR/$AD_EXECUTABLE" ;;
esac

[ -x "$IAD_EXECUTABLE" ] || fail "IAD_EXECUTABLE is not executable: $IAD_EXECUTABLE"

TMPDIR_BASE=${TMPDIR:-/tmp}
TEST_TMP=$(mktemp -d "$TMPDIR_BASE/iad_cli.XXXXXX")
trap 'rm -rf "$TEST_TMP"' EXIT HUP INT TERM

announce() {
    printf '== %s\n' "$*"
}

numeric_actual_file() {
    printf '%s/%s.actual\n' "$TEST_TMP" "$1"
}

run_iad_numeric() {
    label=$1
    shift
    output="$TEST_TMP/$label.out"
    "$IAD_EXECUTABLE" "$@" > "$output" 2>&1 || {
        cat "$output" >&2
        fail "iad failed for $label: $*"
    }
    "$PYTHON" "$CLI_DIR/compare_numeric.py" --extract-last "$output" >> "$ACTUAL_FILE"
}

run_ad_numeric() {
    label=$1
    shift
    output="$TEST_TMP/$label.out"
    "$AD_EXECUTABLE" "$@" > "$output" 2>&1 || {
        cat "$output" >&2
        fail "ad failed for $label: $*"
    }
    "$PYTHON" "$CLI_DIR/compare_numeric.py" --extract-last "$output" >> "$ACTUAL_FILE"
}

compare_numeric_file() {
    expected=$1
    actual=$2
    tolerance=${3:-1e-4}
    "$PYTHON" "$CLI_DIR/compare_numeric.py" --tolerance "$tolerance" "$expected" "$actual"
}

assert_contains() {
    file=$1
    text=$2
    if ! grep -F -- "$text" "$file" >/dev/null 2>&1; then
        echo "Expected to find: $text" >&2
        echo "In file: $file" >&2
        sed -n '1,120p' "$file" >&2
        exit 1
    fi
}

assert_not_exists() {
    path=$1
    [ ! -e "$path" ] || fail "unexpected path exists: $path"
}

assert_exists() {
    path=$1
    [ -e "$path" ] || fail "expected path missing: $path"
}

assert_fails() {
    output="$TEST_TMP/fail.out"
    if "$@" > "$output" 2>&1; then
        cat "$output" >&2
        fail "command unexpectedly succeeded: $*"
    fi
}
