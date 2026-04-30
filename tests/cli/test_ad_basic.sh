#!/bin/sh

. "$(dirname "$0")/lib.sh"

[ -x "$AD_EXECUTABLE" ] || fail "AD_EXECUTABLE is not executable: $AD_EXECUTABLE"

ACTUAL_FILE=$(numeric_actual_file ad_basic)
: > "$ACTUAL_FILE"

announce "ad basic cases"
run_ad_numeric ad_base -m -a 0.9 -b 1 -g 0.8 -q 8
run_ad_numeric ad_oblique -m -a 0.9 -b 1 -g 0.8 -i 45 -q 8
run_ad_numeric ad_slides -m -a 0.9 -b 1 -g 0.8 -n 1.4 -s 1.5 -t 1.5 -q 8
compare_numeric_file "$CLI_DIR/expected/ad_basic.txt" "$ACTUAL_FILE" 1e-4

help="$TEST_TMP/ad_help.out"
"$AD_EXECUTABLE" -h > "$help" 2>&1
assert_contains "$help" "Usage:  ad"
assert_contains "$help" "-a #"
assert_contains "$help" "-m"

version="$TEST_TMP/ad_version.out"
"$AD_EXECUTABLE" -v > "$version" 2>&1
assert_contains "$version" "ad 3"

work="$TEST_TMP/ad_file"
cat > "$work.in" <<EOF
0.9 1 0.8 1.0 1.0 1.0 0 0 8
EOF
"$AD_EXECUTABLE" -m -o "$work.out" "$work.in" > "$TEST_TMP/ad_file_stdout.out" 2>&1
assert_exists "$work.out"
assert_contains "$work.out" "0.042"
