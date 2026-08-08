#!/bin/sh

# The Monte Carlo lost-light calculation must treat the two slides
# separately.  It used to read only the top-slide parameters and apply them
# to both faces, so a slide on the bottom was either duplicated or lost
# entirely.  These are relationships between runs rather than frozen numbers,
# since the Monte Carlo values move with the photon count and the seed.

. "$(dirname "$0")/lib.sh"

announce "monte carlo slide handling"

MC_LOST_EXECUTABLE=${MC_LOST_EXECUTABLE:-"$ROOT_DIR/mc_lost"}
if [ ! -x "$MC_LOST_EXECUTABLE" ]; then
    MC_LOST_EXECUTABLE="$ROOT_DIR/src/mc_lost"
fi
if [ ! -x "$MC_LOST_EXECUTABLE" ]; then
    echo "   (skipped: build it with 'make mc_lost')"
    exit 0
fi

# common sample: weakly absorbing, 2 mm slides, 25.4 mm port
mc_common="-a 0.954 -b 1.111 -g 0 -n 1.3382 -t 2.735 -P 25.4 -B 2 -p 200000 -s 12345"

# field 9 is ur1_lost, 10 is ut1_lost, 11 is uru_lost, 12 is utu_lost
mc_lost_field() {
    field=$1
    shift
    # shellcheck disable=SC2086
    "$MC_LOST_EXECUTABLE" $mc_common "$@" -m 2>/dev/null | tail -1 | awk -v f="$field" '{print $f}'
}

none=$(mc_lost_field 9 -N 1 -T 0 -j 1 -k 0)
top=$(mc_lost_field 9 -N 1.5 -T 2 -j 1 -k 0)
bottom=$(mc_lost_field 9 -N 1 -T 0 -j 1.5 -k 2)
both=$(mc_lost_field 9 -N 1.5 -T 2 -j 1.5 -k 2)

[ -n "$none" ] && [ -n "$both" ] || fail "mc_lost produced no output"

# a slide on either face must change the answer
[ "$bottom" != "$none" ] || \
    fail "bottom slide ignored: loss with a bottom slide equals the bare slab ($bottom)"
[ "$top" != "$both" ] || \
    fail "top slide duplicated: loss with one top slide equals two slides ($top)"

# glass on one face must lose more than bare slab, less than glass on both
ordered=$(awk -v n="$none" -v t="$top" -v b="$bottom" -v x="$both" \
    'BEGIN { print (n < t && t < x && n < b && b < x) ? "yes" : "no" }')
[ "$ordered" = "yes" ] || \
    fail "loss not ordered none < one < both: none=$none top=$top bottom=$bottom both=$both"

# a bottom slide spreads transmitted light more than a top slide does, and
# a top slide spreads reflected light more than a bottom slide does
r_top=$(mc_lost_field 9 -N 1.5 -T 2 -j 1 -k 0)
r_bottom=$(mc_lost_field 9 -N 1 -T 0 -j 1.5 -k 2)
t_top=$(mc_lost_field 10 -N 1.5 -T 2 -j 1 -k 0)
t_bottom=$(mc_lost_field 10 -N 1 -T 0 -j 1.5 -k 2)

asym=$(awk -v rt="$r_top" -v rb="$r_bottom" -v tt="$t_top" -v tb="$t_bottom" \
    'BEGIN { print (rt > rb && tb > tt) ? "yes" : "no" }')
[ "$asym" = "yes" ] || \
    fail "slide asymmetry wrong: ur1_lost top=$r_top bottom=$r_bottom, ut1_lost top=$t_top bottom=$t_bottom"

# With a port wide enough that nothing is lost, the Monte Carlo totals must
# agree with adding-doubling -- including for asymmetric slides, which is the
# case the old code could not represent at all.
mc_vs_ad() {
    out=$("$MC_LOST_EXECUTABLE" -a 0.9 -b 1 -g 0 -n 1.4 -t 1 -P 1000 -B 0 \
        -p 400000 -s 12345 "$@" 2>/dev/null)
    mc=$(echo "$out" | grep 'MC Calc')
    ad=$(echo "$out" | grep 'AD Calc')
    # fields: 1-4 are the MC values, 5-6 are the words "MC Calc",
    # 7-10 are the AD values
    echo "$mc $ad" | awk '{
        ok = 1
        for (i = 1; i <= 4; i++) {
            d = $i - $(i + 6)
            if (d < 0) d = -d
            if (d > 0.005) ok = 0
        }
        print (ok ? "yes" : "no")
    }'
}

for spec in "-N 1 -T 0 -j 1 -k 0" "-N 1.5 -T 1 -j 1.5 -k 1" \
            "-N 1.5 -T 1 -j 1 -k 0" "-N 1 -T 0 -j 1.5 -k 1"; do
    # shellcheck disable=SC2086
    [ "$(mc_vs_ad $spec)" = "yes" ] || \
        fail "monte carlo disagrees with adding-doubling for '$spec'"
done

# -G n and -G f turn the sample over between the two measurements.  The
# reflection losses must be untouched by that and the transmission losses
# must not be, because lost light is not reciprocal even though the total
# transmittance is.
sphere_opts="-r 0.4 -t 0.3 -S 1 -1 '203.2 25.4 12.7 1 0.98'"
lost_column() {
    column=$1
    slides=$2
    "$IAD_EXECUTABLE" -r 0.4 -t 0.3 -S 1 -1 "203.2 25.4 12.7 1 0.98" \
        -n 1.4 -N 1.5 -d 1 -D 2 -B 2 -p 2000000 -G "$slides" -x 8 2>&1 >/dev/null |
        grep -E '^[[:space:]]+1[[:space:]]' | tail -1 | awk -v c="$column" '{print $c}'
}

for pair in "t n" "b f"; do
    # shellcheck disable=SC2086
    set -- $pair
    plain=$1
    flipped=$2

    r_plain=$(lost_column 12 "$plain")
    r_flipped=$(lost_column 12 "$flipped")
    t_plain=$(lost_column 14 "$plain")
    t_flipped=$(lost_column 14 "$flipped")

    [ "$r_plain" = "$r_flipped" ] || \
        fail "-G $flipped changed the reflection loss from -G $plain: $r_plain vs $r_flipped"
    [ "$t_plain" != "$t_flipped" ] || \
        fail "-G $flipped did not flip the sample for transmission: both give $t_plain"
done
