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
# -w and -W override the wall reflectance that an .rxt file already gave, so
# they are only meaningful with a file.  Combined with -1 or -2 they have
# nothing to override and must be refused.
cp "$ROOT_DIR/tests/rxt/1_sphere/vio_A.rxt" "$TEST_TMP/wall.rxt"
wall_musp() {
    "$IAD_EXECUTABLE" -M 0 -l '650 660' "$@" -o "$TEST_TMP/wall.txt" \
        "$TEST_TMP/wall.rxt" >/dev/null 2>&1 || fail "iad failed for -w test: $*"
    grep -v '^#' "$TEST_TMP/wall.txt" | head -1 | awk '{print $7}'
}
wall_default=$(wall_musp)
wall_low=$(wall_musp -w 0.80)
wall_high=$(wall_musp -w 0.99)
[ "$wall_default" != "$wall_low" ] && [ "$wall_low" != "$wall_high" ] || \
    fail "-w did not override the wall reflectance in the file: $wall_default $wall_low $wall_high"

assert_fails "$IAD_EXECUTABLE" -r 0.2 -t 0.01 -S 1 -1 "100 15 13 2 0.95" -w 0.95
assert_fails "$IAD_EXECUTABLE" -r 0.2 -t 0.01 -S 2 -1 "100 15 13 2 0.95" -2 "100 15 0 2 0.95" -W 0.95
assert_fails "$IAD_EXECUTABLE" -r 0.2 -t 0.01 -S 2 -1 "100 15 13 2 0.95" -2 "100 15 0 2 0.95" -w 0.95
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

# a failed lost-light correction is 'm' and names the Monte Carlo loop.  This
# needs a real sphere: it used to rely on the default geometry, whose sample
# port is zero, which is now refused outright.  The geometry below is combo_0's
# -- a 44.45 mm sample port with a 6.5 mm beam through a 6.5 mm entrance port.
mc_stalled="$TEST_TMP/iad_mc_stalled.out"
"$IAD_EXECUTABLE" -r 0.72 -t 0.26 -n 1.365 -N 1.5875 -d 2 -D 1.075 -B 6.5 -S 1 \
    -1 '203.2 44.45 6.5 1 0.97' -2 '203.2 44.45 0 1 0.97' > "$mc_stalled" 2>&1
assert_matches "$mc_stalled" "0\.7200.*0\.2600.*m"
assert_contains "$mc_stalled" "lost-light correction did not settle"

# A slide removed by -G must not keep the absorption set by -E.  These are
# relationships between runs rather than frozen numbers, so they stay valid
# if the underlying values ever shift.
slab_opts="-a 0.9 -b 1 -g 0 -n 1.4 -N 1.5 -d 1 -D 1 -q 32"
t_total() {
    "$IAD_EXECUTABLE" $slab_opts "$@" -z 2>/dev/null | grep 'T total' | awk '{print $NF}'
}

no_slide_e0=$(t_total -G 0 -E 0)
no_slide_e5=$(t_total -G 0 -E 0.5)
[ "$no_slide_e0" = "$no_slide_e5" ] || \
    fail "-G 0 says no slides, but -E 0.5 still absorbed: $no_slide_e0 vs $no_slide_e5"

one_slide_e5=$(t_total -G t -E 0.5)
two_slide_e5=$(t_total -G 2 -E 0.5)
[ "$one_slide_e5" != "$two_slide_e5" ] || \
    fail "-G t -E 0.5 matches -G 2 -E 0.5; absorption applied to an absent slide"

# ...but an index-matched absorbing slide is legitimate and must survive
matched=$(t_total -E 0.5)
[ "$matched" = "$two_slide_e5" ] || \
    fail "default two-slide absorption changed: $matched vs $two_slide_e5"

# The beam has to fit through the ports it passes through.
beam_status() {
    "$IAD_EXECUTABLE" -r 0.3 -t 0.1 -S 1 "$@" 2>/dev/null |
        grep -v '^#' | tail -1 | awk '{print $NF}'
}
[ "$(beam_status -1 '203.2 44.45 6.3 1 0.97' -B 6.5)" = "P" ] || \
    fail "beam wider than the entrance port was not rejected"
[ "$(beam_status -1 '203.2 5 12.7 1 0.97' -B 6.5)" = "P" ] || \
    fail "beam wider than the sample port was not rejected"
[ "$(beam_status -1 '203.2 44.45 12.7 1 0.97' -B 6.5)" = "*" ] || \
    fail "a beam that fits through both ports was rejected"

# With two spheres the sample sits in one hole, so both blocks must give it
# the same diameter.  The port areas are fractions of their own sphere, so
# this has to compare diameters -- the second case below has area fractions
# differing by a factor of four yet describes the same physical port.
ports_status() {
    "$IAD_EXECUTABLE" -r 0.4 -t 0.18 -S 2 "$@" 2>/dev/null |
        grep -v '^#' | tail -1 | awk '{print $NF}'
}
ports_status_n() {
    n=$1
    shift
    "$IAD_EXECUTABLE" -r 0.4 -t 0.18 -S "$n" "$@" 2>/dev/null |
        grep -v '^#' | tail -1 | awk '{print $NF}'
}
mismatch="$TEST_TMP/ports_mismatch.out"
"$IAD_EXECUTABLE" -r 0.4 -t 0.18 -S 2 \
    -1 '203.2 31.75 6.35 3.18 0.975' -2 '203.2 15.9 6.35 3.18 0.975' \
    > "$mismatch" 2>&1
assert_contains "$mismatch" "the two spheres disagree about the sample port"

[ "$(ports_status -1 '203.2 31.75 6.35 3.18 0.975' -2 '203.2 15.9 6.35 3.18 0.975')" = "P" ] || \
    fail "two spheres with different sample ports were accepted"
[ "$(ports_status -1 '203.2 31.75 6.35 3.18 0.975' -2 '101.6 31.75 6.35 3.18 0.975')" = "*" ] || \
    fail "same sample port in differently sized spheres was rejected"

# "One sphere" means one sphere at a time, not one sphere in the experiment:
# the reflectance may be measured on one sphere and the transmittance on
# another of a different size with a different sample port.  That must be
# allowed.  The very same geometry has to be refused when both spheres are
# declared present at once, because then the sample faces both through a
# single hole.
[ "$(ports_status_n 1 -1 '203.2 25.4 12.7 1 0.97' -2 '101.6 12.7 0 1 0.97')" = "*" ] || \
    fail "one sphere at a time: different spheres for R and T were rejected"
[ "$(ports_status_n 2 -1 '203.2 25.4 12.7 1 0.97' -2 '101.6 12.7 0 1 0.97')" = "P" ] || \
    fail "two spheres at once: mismatched sample ports were accepted"

# The transmission sphere's exit port may be any size including none.  Closed
# over, the direct beam lands on the sphere wall.  Open but narrower than the
# beam cannot work, since only part of the direct beam would reach the
# standard sitting in the port.
exit_status() {
    port=$1
    "$IAD_EXECUTABLE" -r 0.3 -t 0.1 -S 1 -1 '203.2 25.4 12.7 1 0.97' \
        -2 "203.2 25.4 $port 1 0.97" -B 6.5 2>/dev/null |
        grep -v '^#' | tail -1 | awk '{print $NF}'
}
[ "$(exit_status 0)" = "*" ] || \
    fail "a closed exit port was rejected; zero must mean the beam hits the wall"
[ "$(exit_status 3)" = "P" ] || \
    fail "an exit port narrower than the beam was accepted"
[ "$(exit_status 12.7)" = "*" ] || \
    fail "an exit port wider than the beam was rejected"

narrow_exit="$TEST_TMP/exit_narrow.out"
"$IAD_EXECUTABLE" -r 0.3 -t 0.1 -S 1 -1 '203.2 25.4 12.7 1 0.97' \
    -2 '203.2 25.4 3 1 0.97' -B 6.5 > "$narrow_exit" 2>&1
assert_contains "$narrow_exit" "beam is wider than the exit port"

assert_fails "$IAD_EXECUTABLE" -a 2 -r 0.2
assert_fails "$IAD_EXECUTABLE" -q 5 -r 0.2
assert_fails "$IAD_EXECUTABLE" -S 3 -r 0.2
assert_fails "$IAD_EXECUTABLE" "$TEST_TMP/missing.rxt"
