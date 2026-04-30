#!/bin/sh

. "$(dirname "$0")/lib.sh"

ACTUAL_FILE=$(numeric_actual_file iad_two_sphere_no_mc)
: > "$ACTUAL_FILE"

announce "iad two-sphere deterministic cases"
run_iad_numeric two_reflectance -V 0 -r 0.4 -S 2 -M 0 -1 "200 13 13 2 0.95" -2 "200 13 0 2 0.95"
run_iad_numeric two_rt -V 0 -r 0.4 -t 0.1 -S 2 -M 0 -1 "200 13 13 2 0.95" -2 "200 13 0 2 0.95"
run_iad_numeric two_rtu -V 0 -r 0.4 -t 0.1 -u 0.002 -S 2 -M 0 -1 "200 13 13 2 0.95" -2 "200 13 0 2 0.95"
run_iad_numeric two_lower_rt -V 0 -r 0.2 -t 0.1 -u 0.0049787 -S 2 -M 0 -1 "200 13 13 2 0.95" -2 "200 13 0 2 0.95"

compare_numeric_file "$CLI_DIR/expected/iad_two_sphere_no_mc.txt" "$ACTUAL_FILE" 8e-4
