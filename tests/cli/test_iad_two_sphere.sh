#!/bin/sh

. "$(dirname "$0")/lib.sh"

announce "iad two-sphere MC smoke case"
out="$TEST_TMP/two_sphere_mc.out"
"$IAD_EXECUTABLE" -V 0 -r 0.2 -t 0.1 -u 0.0049787 -S 2 -M 1 -p 1000 \
    -1 "200 13 13 2 0.95" -2 "200 13 0 2 0.95" > "$out" 2>&1
"$PYTHON" "$CLI_DIR/compare_numeric.py" --extract-last "$out" > "$TEST_TMP/two_sphere_mc.numeric"
assert_contains "$TEST_TMP/two_sphere_mc.numeric" "0.2"
