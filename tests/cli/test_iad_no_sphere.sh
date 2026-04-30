#!/bin/sh

. "$(dirname "$0")/lib.sh"

announce "iad no-sphere fixture cases"

cp "$ROOT_DIR/tests/rxt/0_sphere/basic_A.rxt" "$TEST_TMP/basic_A.rxt"
cp "$ROOT_DIR/tests/rxt/0_sphere/sample_A.rxt" "$TEST_TMP/sample_A.rxt"

(
    cd "$TEST_TMP"
    "$IAD_EXECUTABLE" -M 0 -q 4 basic_A.rxt > basic_A.stdout 2>&1
    "$IAD_EXECUTABLE" -M 0 -q 4 -o sample_A.out sample_A.rxt > sample_A.stdout 2>&1
)

assert_exists "$TEST_TMP/basic_A.txt"
assert_exists "$TEST_TMP/sample_A.out"
assert_contains "$TEST_TMP/basic_A.txt" "# Inverse Adding-Doubling"
assert_contains "$TEST_TMP/sample_A.out" "*"
