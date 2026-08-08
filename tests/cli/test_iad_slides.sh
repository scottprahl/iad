#!/bin/sh

# Round trip through every slide configuration.  For each one, -z computes
# R, T and unscattered T from known optical properties, and those three
# numbers are fed straight back into the inverse.  The inverse must recover
# what we started with.
#
# Sample: a=0.8, b=5, g=0.9, index 1.5, 1 mm thick, so
#     mu_a   = (1-a) b/d = 1.0  1/mm
#     mu_s'  = a b/d (1-g) = 0.4  1/mm
#     g      = 0.9
# Slides: index 1.4, 1 mm thick, absorption from -E.
#
# No spheres are used, so no Monte Carlo lost-light estimate enters and every
# number is deterministic.  -G n and -G f are omitted because they are
# identical to -G t and -G b for reflectance and transmittance; the flip only
# shows up in lost light, which test_mc_slides.sh covers.

. "$(dirname "$0")/lib.sh"

announce "slide configurations (forward then inverse)"

geom="-n 1.5 -N 1.4 -d 1 -D 1 -q 8"

want_mu_a=1.0
want_mu_sp=0.4
want_g=0.9
want_b=5.0
thickness=1.0
tolerance=0.02

round_trip() {
    slides=$1
    absorption=$2

    forward="$TEST_TMP/fwd_${slides}_${absorption}.out"
    # shellcheck disable=SC2086
    "$IAD_EXECUTABLE" -z -a 0.8 -b 5 -g 0.9 $geom -G "$slides" -E "$absorption" \
        > "$forward" 2>&1 || {
        cat "$forward" >&2
        fail "forward calculation failed for -G $slides -E $absorption"
    }

    R=$(grep 'R total'       "$forward" | awk '{print $NF}')
    T=$(grep 'T total'       "$forward" | awk '{print $NF}')
    Tu=$(grep 'T unscattered' "$forward" | awk '{print $NF}')

    [ -n "$R" ] && [ -n "$T" ] && [ -n "$Tu" ] || \
        fail "could not read R, T and Tu for -G $slides -E $absorption"

    inverse="$TEST_TMP/inv_${slides}_${absorption}.out"
    # shellcheck disable=SC2086
    "$IAD_EXECUTABLE" -r "$R" -t "$T" -u "$Tu" $geom -G "$slides" -E "$absorption" \
        > "$inverse" 2>/dev/null
    grep -v '^#' "$inverse" | tail -1
}

check_recovered() {
    slides=$1
    absorption=$2

    row=$(round_trip "$slides" "$absorption")
    status=$(echo "$row" | awk '{print $NF}')
    [ "$status" = "*" ] || \
        fail "-G $slides -E $absorption did not converge (status '$status'): $row"

    # iad reports mu_a, mu_s' and g.  The optical thickness that the forward
    # calculation used has to come back too, so rebuild it:
    #     mu_s = mu_s' / (1 - g)      b = (mu_a + mu_s) * d
    echo "$row" | awk -v a="$want_mu_a" -v s="$want_mu_sp" -v g="$want_g" \
        -v b="$want_b" -v d="$thickness" \
        -v tol="$tolerance" -v label="-G $slides -E $absorption" '
    function off(got, want) { return (want == 0) ? got : (got - want) / want }
    function abs(x) { return x < 0 ? -x : x }
    {
        mu_a = $6; mu_sp = $7; gg = $8
        if (gg >= 1.0) {
            printf "FAIL: %s returned g=%s, cannot recover b\n", label, gg > "/dev/stderr"
            exit 1
        }
        got_b = (mu_a + mu_sp / (1.0 - gg)) * d

        if (abs(off(mu_a, a)) > tol || abs(off(mu_sp, s)) > tol ||
            abs(off(gg, g)) > tol || abs(off(got_b, b)) > tol) {
            printf "FAIL: %s recovered mu_a=%s mu_s'\''=%s g=%s b=%.4f, wanted %s %s %s %s\n",
                   label, mu_a, mu_sp, gg, got_b, a, s, g, b > "/dev/stderr"
            exit 1
        }
    }' || exit 1
}

check_recovered 0 0
check_recovered 2 0
check_recovered t 0
check_recovered b 0
check_recovered 0 0.5
check_recovered 2 0.5
check_recovered t 0.5
check_recovered b 0.5
