#!/bin/sh

. "$(dirname "$0")/lib.sh"

mode=${1:-full}
ACTUAL_FILE=$(numeric_actual_file iad_basic)
: > "$ACTUAL_FILE"

announce "iad basic numeric cases ($mode)"

run_common_cases() {
    run_iad_numeric r0 -V 0 -r 0
    run_iad_numeric r1 -V 0 -r 1
    run_iad_numeric r04 -V 0 -r 0.4
    run_iad_numeric r04_t01 -V 0 -r 0.4 -t 0.1
    run_iad_numeric r04_t01_u0002 -V 0 -r 0.4 -t 0.1 -u 0.002
    run_iad_numeric one_sphere_minimal -V 0 -r 0.4 -t 0.01 -d 1 -M 0 -S 1 -1 "200 13 13 2 0.95"
}

run_common_cases

if [ "$mode" = "--veryshort" ] || [ "$mode" = "veryshort" ]; then
    compare_numeric_file "$CLI_DIR/expected/iad_basic_veryshort.txt" "$ACTUAL_FILE" 8e-4
    exit 0
fi

run_iad_numeric r04_t01_u0049787 -V 0 -r 0.4 -t 0.1 -u 0.049787

run_iad_numeric n15_r04 -V 0 -r 0.4 -n 1.5
run_iad_numeric n15_r04_t01 -V 0 -r 0.4 -t 0.1 -n 1.5
run_iad_numeric n15_r04_t01_u0002 -V 0 -r 0.4 -t 0.1 -u 0.002 -n 1.5
run_iad_numeric n15_r04_t01_u0045884 -V 0 -r 0.4 -t 0.1 -u 0.045884 -n 1.5

run_iad_numeric slides_r04 -V 0 -r 0.4 -n 1.4 -N 1.5
run_iad_numeric slides_r04_t01 -V 0 -r 0.4 -t 0.1 -n 1.4 -N 1.5
run_iad_numeric slides_r04_t01_u0002 -V 0 -r 0.4 -t 0.1 -u 0.002 -n 1.4 -N 1.5
run_iad_numeric slides_r04_t01_u0045884 -V 0 -r 0.4 -t 0.1 -u 0.045884 -n 1.4 -N 1.5

run_iad_numeric top_slide_r04 -V 0 -r 0.4 -n 1.4 -N 1.5 -G t
run_iad_numeric top_slide_r04_t01 -V 0 -r 0.4 -t 0.1 -n 1.4 -N 1.5 -G t
run_iad_numeric top_slide_r04_t01_u0002 -V 0 -r 0.4 -t 0.1 -u 0.002 -n 1.4 -N 1.5 -G t
run_iad_numeric top_slide_r04_t01_u0045884 -V 0 -r 0.4 -t 0.1 -u 0.045884 -n 1.4 -N 1.5 -G t

run_iad_numeric bottom_slide_r04 -V 0 -r 0.4 -n 1.4 -N 1.5 -G b
run_iad_numeric bottom_slide_r04_t01 -V 0 -r 0.4 -t 0.1 -n 1.4 -N 1.5 -G b
run_iad_numeric bottom_slide_r04_t01_u0002 -V 0 -r 0.4 -t 0.1 -u 0.002 -n 1.4 -N 1.5 -G b
run_iad_numeric bottom_slide_r04_t01_u0045884 -V 0 -r 0.4 -t 0.1 -u 0.045884 -n 1.4 -N 1.5 -G b

run_iad_numeric absorbing_slide_1 -V 0 -r 0.0000000 -t 0.135335 -E 0.5
run_iad_numeric absorbing_slide_2 -V 0 -r 0.0249268 -t 0.155858 -E 0.5
run_iad_numeric absorbing_slide_3 -V 0 -r 0.0520462 -t 0.134587 -E 0.5 -n 1.5 -N 1.5

run_iad_numeric g_r04 -V 0 -r 0.4 -g 0.9
run_iad_numeric g_r04_t01 -V 0 -r 0.4 -t 0.1 -g 0.9
run_iad_numeric g_n15_r04 -V 0 -r 0.4 -g 0.9 -n 1.5
run_iad_numeric g_n15_r04_t01 -V 0 -r 0.4 -t 0.1 -g 0.9 -n 1.5
run_iad_numeric g_slides_r04 -V 0 -r 0.4 -g 0.9 -n 1.4 -N 1.5
run_iad_numeric g_slides_r04_t01 -V 0 -r 0.4 -t 0.1 -g 0.9 -n 1.4 -N 1.5

run_iad_numeric a09_r04 -V 0 -r 0.4 -a 0.9
run_iad_numeric a09_r04_t01 -V 0 -r 0.4 -t 0.1 -a 0.9
run_iad_numeric a09_n15_r04_t01 -V 0 -r 0.4 -t 0.1 -a 0.9 -n 1.5
run_iad_numeric a095_slides_r04_t01 -V 0 -r 0.4 -t 0.1 -a 0.95 -n 1.4 -N 1.5

run_iad_numeric b3_r04_t01 -V 0 -r 0.4 -t 0.1 -b 3
run_iad_numeric b3_n15_r04_t01 -V 0 -r 0.4 -t 0.1 -b 3 -n 1.5
run_iad_numeric b3_slides_r04_t01 -V 0 -r 0.4 -t 0.1 -b 3 -n 1.4 -N 1.5

run_iad_numeric F30_r04 -V 0 -r 0.4 -F 30
run_iad_numeric F30_r04_t01 -V 0 -r 0.4 -t 0.1 -F 30
run_iad_numeric F30_n15_r04 -V 0 -r 0.4 -F 30 -n 1.5
run_iad_numeric F30_n15_r04_t01 -V 0 -r 0.4 -t 0.1 -F 30 -n 1.5
run_iad_numeric F30_slides_r04 -V 0 -r 0.4 -F 30 -n 1.4 -N 1.5
run_iad_numeric F30_slides_r04_t01 -V 0 -r 0.4 -t 0.1 -F 30 -n 1.4 -N 1.5

run_iad_numeric A06_r03 -V 0 -r 0.3 -A 0.6
run_iad_numeric A06_r03_t01 -V 0 -r 0.3 -t 0.1 -A 0.6
run_iad_numeric A06_n15_r03 -V 0 -r 0.3 -A 0.6 -n 1.5
run_iad_numeric A06_n15_r03_t01 -V 0 -r 0.3 -t 0.1 -A 0.6 -n 1.5
run_iad_numeric A06_slides_r03 -V 0 -r 0.3 -A 0.6 -n 1.4 -N 1.5
run_iad_numeric A06_slides_r03_t01 -V 0 -r 0.3 -t 0.1 -A 0.6 -n 1.4 -N 1.5

run_iad_numeric A06_g06_r03 -V 0 -r 0.3 -A 0.6 -g 0.6
run_iad_numeric A06_g06_n15_r03 -V 0 -r 0.3 -A 0.6 -g 0.6 -n 1.5
run_iad_numeric A06_g06_slides_r03 -V 0 -r 0.3 -A 0.6 -g 0.6 -n 1.4 -N 1.5

run_iad_numeric F2_g05_r03 -V 0 -r 0.3 -F 2.0 -g 0.5
run_iad_numeric F2_g05_n15_r03 -V 0 -r 0.3 -F 2.0 -g 0.5 -n 1.5
run_iad_numeric F2_g05_slides_r03 -V 0 -r 0.3 -F 2.0 -g 0.5 -n 1.4 -N 1.5

run_iad_numeric oblique_absorbing -V 0 -i 60 -r 0.00000 -t 0.13691
run_iad_numeric oblique_mixed -V 0 -i 60 -r 0.14932 -t 0.23181
run_iad_numeric oblique_scattering -V 0 -i 60 -r 0.61996 -t 0.30605

compare_numeric_file "$CLI_DIR/expected/iad_basic.txt" "$ACTUAL_FILE" 2e-3
