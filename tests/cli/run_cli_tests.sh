#!/bin/sh

set -eu

case "${1:-full}" in
    veryshort)
        "$(dirname "$0")/test_iad_basic.sh" --veryshort
        ;;
    basic)
        "$(dirname "$0")/test_iad_basic.sh"
        "$(dirname "$0")/test_iad_forward.sh"
        "$(dirname "$0")/test_iad_no_sphere.sh"
        "$(dirname "$0")/test_iad_one_sphere_no_mc.sh"
        "$(dirname "$0")/test_iad_two_sphere_no_mc.sh"
        "$(dirname "$0")/test_iad_files.sh"
        "$(dirname "$0")/test_iad_options.sh"
        "$(dirname "$0")/test_ad_basic.sh"
        ;;
    full)
        "$(dirname "$0")/run_cli_tests.sh" basic
        "$(dirname "$0")/test_iad_one_sphere.sh"
        "$(dirname "$0")/test_iad_two_sphere.sh"
        ;;
    batch)
        "$(dirname "$0")/test_iad_files.sh" --batch
        ;;
    *)
        echo "usage: $0 [veryshort|basic|full|batch]" >&2
        exit 2
        ;;
esac
