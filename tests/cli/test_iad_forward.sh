#!/bin/sh

. "$(dirname "$0")/lib.sh"

announce "iad forward calculations"

out="$TEST_TMP/forward.out"
"$IAD_EXECUTABLE" -z -a 0.9 -b 1 -g 0.8 -q 8 > "$out" 2>&1
assert_contains "$out" "R total         = 0.042"
assert_contains "$out" "T total         = 0.846"
assert_contains "$out" "T unscattered   = 0.368"

"$IAD_EXECUTABLE" -z -a 0.9 -b 1 -g 0.8 -i 45 -q 12 > "$out" 2>&1
assert_contains "$out" "cos(theta incident) = 0.707"
assert_contains "$out" "Calculated quantities"

"$IAD_EXECUTABLE" -z -a 0.95 -b 2 -g 0.9 -n 1.4 -N 1.5 -G 2 -q 8 > "$out" 2>&1
assert_contains "$out" "sample index        = 1.400"
assert_contains "$out" "top slide index     = 1.500"
assert_contains "$out" "bottom slide index  = 1.500"

"$IAD_EXECUTABLE" -z -V 0 -a 0.9 -b 1 -g 0.8 -q 8 > "$out" 2>&1
[ ! -s "$out" ] || fail "iad -z -V 0 should not write stdout"
