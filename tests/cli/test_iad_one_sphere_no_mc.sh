#!/bin/sh

. "$(dirname "$0")/lib.sh"

ACTUAL_FILE=$(numeric_actual_file iad_one_sphere_no_mc)
: > "$ACTUAL_FILE"

announce "iad one-sphere deterministic cases"
run_iad_numeric one_reflectance -V 0 -r 0.4 -S 1 -M 0 -1 "200 13 13 2 0.95"
run_iad_numeric one_rt -V 0 -r 0.4 -t 0.1 -S 1 -M 0 -1 "200 13 13 2 0.95"
run_iad_numeric one_rtu -V 0 -r 0.4 -t 0.1 -u 0.002 -S 1 -M 0 -1 "200 13 13 2 0.95"
run_iad_numeric one_baffles -V 0 -r 0.2 -t 0.2 -u 0.0049787 -S 1 -M 0 -1 "100 13 13 2 0.95" -2 "100 13 0 2 0.95" -H 3

compare_numeric_file "$CLI_DIR/expected/iad_one_sphere_no_mc.txt" "$ACTUAL_FILE" 2e-3
