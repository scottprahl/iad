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
assert_matches "$version" "^iad [0-9]+-[0-9]+-[0-9]+ "

run_iad_numeric options_refs -V 0 -r 0.2 -t 0.01 -M 0 -S 1 \
    -1 "100 15 13 2 0.95" -2 "100 15 13 2 0.95" -T 0.5 -R 0.5
run_iad_numeric options_fractions -V 0 -r 0.4 -t 0.1 -u 0.002 -M 0 -c 0.2 -C 0.3 -f 0.5
run_iad_numeric options_constraints -V 0 -r 0.3 -t 0.1 -M 0 -A 0.6 -g 0.6 -F "P 500 2.0 -1.0" -L 500 -d 1 -D 0 -B 10
run_iad_numeric options_wall_defaults -V 0 -r 0.2 -t 0.01 -M 0 -S 1 -w 0.95 -W 0.95
run_iad_numeric options_search_debug -V 0 -r 0.2 -t 0.01 -M 0 -S 1 -s 2 -H 3 -x 0 \
    -1 "100 15 13 2 0.95" -2 "100 15 13 2 0.95"

agrid_bag="$TEST_TMP/agrid_bag.out"
"$IAD_EXECUTABLE" -V 0 -x 2 -r 0.4 -t 0.1 -F 30 > "$agrid_bag" 2>&1
assert_contains "$agrid_bag" "AGRID: Filling BaG grid"

agrid_bsg="$TEST_TMP/agrid_bsg.out"
"$IAD_EXECUTABLE" -V 0 -x 2 -r 0.3 -t 0.1 -A 0.6 > "$agrid_bsg" 2>&1
assert_contains "$agrid_bsg" "AGRID: Filling BsG grid"

cp "$ROOT_DIR/tests/rxt/1_sphere/vio_A.rxt" "$TEST_TMP/vio_A.rxt"
(
    cd "$TEST_TMP"
    "$IAD_EXECUTABLE" -M 0 -q 4 -J -l "650 700" vio_A.rxt > grid.out 2>&1
)
assert_exists "$TEST_TMP/vio_A.txt"

# M_R + M_T > 1 violates energy conservation when no spheres are used, so the
# row is flagged '!' rather than silently reporting a bogus fit.
too_much="$TEST_TMP/iad_too_much_light.out"
"$IAD_EXECUTABLE" -r 0.4 -t 0.7 > "$too_much" 2>&1
assert_matches "$too_much" "0\.4000.*0\.7000.*!"
# a single command-line measurement names its one error rather than
# printing the whole legend
assert_contains "$too_much" "Failed Search, M_R + M_T exceeds 1"

# a sum of exactly 1 is the lossless limit and must still invert
lossless="$TEST_TMP/iad_lossless.out"
"$IAD_EXECUTABLE" -r 0.4 -t 0.6 > "$lossless" 2>&1
assert_matches "$lossless" "0\.4000.*0\.6000.*\*"

# a search that stops early is 'x', not the old misleading "too many iterations"
stalled="$TEST_TMP/iad_stalled.out"
"$IAD_EXECUTABLE" -r 0.9 -t 0.05 -n 1.5 > "$stalled" 2>&1
assert_matches "$stalled" "0\.9000.*0\.0500.*x"
assert_contains "$stalled" "stopped early without matching the measurements"

# a failed lost-light correction is 'm' and names the Monte Carlo loop
mc_stalled="$TEST_TMP/iad_mc_stalled.out"
"$IAD_EXECUTABLE" -r 0.7 -t 0.28 -n 1.33 -S 1 > "$mc_stalled" 2>&1
assert_matches "$mc_stalled" "0\.7000.*0\.2800.*m"
assert_contains "$mc_stalled" "lost-light correction did not settle"

assert_fails "$IAD_EXECUTABLE" -a 2 -r 0.2
assert_fails "$IAD_EXECUTABLE" -q 5 -r 0.2
assert_fails "$IAD_EXECUTABLE" -S 3 -r 0.2
assert_fails "$IAD_EXECUTABLE" "$TEST_TMP/missing.rxt"
