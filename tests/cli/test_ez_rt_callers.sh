#!/bin/sh

# The ez_RT family (ez_RT, ez_RT_unscattered, ez_RT_Cone, ez_RT_Oblique) takes
# indices of refraction but no slide optical depths, and each one sets both to
# zero internally.  The slides refract but never absorb.  That is deliberate --
# they are small interfaces for callers that cannot pass a structure -- but it
# means their answers are wrong for a sample with absorbing slides, and the
# signature gives no hint of it.
#
# So nothing may call them from anywhere that slide absorption can reach.  Two
# callers are allowed:
#
#   ad_layers_test  -- a standalone test that uses clear slides throughout
#   mc_lost_main    -- its slide absorption is fixed at zero and cannot be set
#                      from the mc_lost command line
#
# A new caller elsewhere -- iad_main, iad_calc, iad_pub, ad_main -- would
# silently drop slide absorption, exactly the bug that made iad -z disagree
# with the inverse.  This test fails if one appears.

. "$(dirname "$0")/lib.sh"

announce "ez_RT callers cannot see absorbing slides"

SRC_DIR="$ROOT_DIR/src"
[ -d "$SRC_DIR" ] || fail "cannot find src directory at $SRC_DIR"

# where the family is defined, and who is allowed to call it
defines='ad_prime|ad_cone'
allowed='ad_layers_test|mc_lost_main'

offenders="$TEST_TMP/ez_rt_offenders"
: > "$offenders"

for f in "$SRC_DIR"/*.w "$SRC_DIR"/*.c; do
    [ -f "$f" ] || continue
    base=$(basename "$f")
    stem=${base%.*}

    echo "$stem" | grep -Eq "^($defines)$" && continue
    echo "$stem" | grep -Eq "^($allowed)$" && continue

    # a call is the name followed by '(' -- skip prototypes and header text
    if grep -nE 'ez_RT[A-Za-z_]*[[:space:]]*\(' "$f" | grep -vE 'Prototype|@<|^\s*\*' > "$TEST_TMP/hits" 2>/dev/null; then
        if [ -s "$TEST_TMP/hits" ]; then
            sed "s|^|$base:|" "$TEST_TMP/hits" >> "$offenders"
        fi
    fi
done

if [ -s "$offenders" ]; then
    echo "The ez_RT family ignores slide absorption; these callers are not allowed:" >&2
    cat "$offenders" >&2
    echo "Use RT with a filled-in AD_slab_type, or Sp_mu_RT for the unscattered part." >&2
    fail "ez_RT called from a file where slide absorption can be set"
fi

# The allowance for mc_lost_main rests on its slide absorption being fixed at
# zero.  If a command-line option for it ever appears, the comparison against
# adding-doubling there becomes wrong and this allowance must be revisited.
if grep -q "double b_slide = 0;" "$SRC_DIR/mc_lost_main.c"; then
    :
else
    fail "mc_lost_main.c no longer fixes b_slide at zero; recheck its ez_RT_Oblique call"
fi
