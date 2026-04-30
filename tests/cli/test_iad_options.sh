#!/bin/sh

. "$(dirname "$0")/lib.sh"

announce "iad option parsing and behavior"
ACTUAL_FILE=$(numeric_actual_file iad_options)
: > "$ACTUAL_FILE"

help="$TEST_TMP/iad_help.out"
"$IAD_EXECUTABLE" -h > "$help" 2>&1
assert_contains "$help" "Usage:  iad"
assert_contains "$help" "-1 '# # # # #'"
assert_contains "$help" "-z"

version="$TEST_TMP/iad_version.out"
"$IAD_EXECUTABLE" -v > "$version" 2>&1
assert_contains "$version" "iad 3"

run_iad_numeric options_refs -V 0 -r 0.2 -t 0.01 -M 0 -S 1 \
    -1 "100 15 13 2 0.95" -2 "100 15 13 2 0.95" -T 0.5 -R 0.5
run_iad_numeric options_fractions -V 0 -r 0.4 -t 0.1 -u 0.002 -M 0 -c 0.2 -C 0.3 -f 0.5
run_iad_numeric options_constraints -V 0 -r 0.3 -t 0.1 -M 0 -A 0.6 -g 0.6 -F "P 500 2.0 -1.0" -L 500 -d 1 -D 0 -B 10
run_iad_numeric options_wall_defaults -V 0 -r 0.2 -t 0.01 -M 0 -S 1 -w 0.95 -W 0.95
run_iad_numeric options_search_debug -V 0 -r 0.2 -t 0.01 -M 0 -S 1 -s 2 -H 3 -x 0 \
    -1 "100 15 13 2 0.95" -2 "100 15 13 2 0.95"

cp "$ROOT_DIR/tests/rxt/1_sphere/vio_A.rxt" "$TEST_TMP/vio_A.rxt"
(
    cd "$TEST_TMP"
    "$IAD_EXECUTABLE" -M 0 -q 4 -J -l "650 700" vio_A.rxt > grid.out 2>&1
)
assert_exists "$TEST_TMP/vio_A.txt"

assert_fails "$IAD_EXECUTABLE" -a 2 -r 0.2
assert_fails "$IAD_EXECUTABLE" -q 5 -r 0.2
assert_fails "$IAD_EXECUTABLE" -S 3 -r 0.2
assert_fails "$IAD_EXECUTABLE" "$TEST_TMP/missing.rxt"
