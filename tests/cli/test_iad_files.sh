#!/bin/sh

. "$(dirname "$0")/lib.sh"

announce "iad file workflow cases"

copy_fixture() {
    src=$1
    dest=$2
    cp "$ROOT_DIR/tests/rxt/$src" "$TEST_TMP/$dest"
}

run_file_case() {
    input=$1
    output=$2
    shift 2
    "$IAD_EXECUTABLE" -M 0 -q 4 "$@" "$input" > "$TEST_TMP/file_case.out" 2>&1 || {
        cat "$TEST_TMP/file_case.out" >&2
        fail "iad file case failed for $input"
    }
    assert_exists "$output"
    assert_contains "$output" "*"
}

if [ "${1:-}" = "--batch" ]; then
    copy_fixture "0_sphere/basic_A.rxt" "basic_A.rxt"
    copy_fixture "1_sphere/vio_A.rxt" "vio_A.rxt"
    copy_fixture "2_sphere/sample_E.rxt" "sample_E.rxt"
    (
        cd "$TEST_TMP"
        "$IAD_EXECUTABLE" -M 0 -q 4 basic_A.rxt > batch_0.out 2>&1
        "$IAD_EXECUTABLE" -M 0 -q 4 vio_A.rxt > batch_1.out 2>&1
        "$IAD_EXECUTABLE" -M 0 -q 4 sample_E.rxt > batch_2.out 2>&1
    )
    assert_exists "$TEST_TMP/basic_A.txt"
    assert_exists "$TEST_TMP/vio_A.txt"
    assert_exists "$TEST_TMP/sample_E.txt"
    exit 0
fi

copy_fixture "0_sphere/basic_A.rxt" "basic_A.rxt"
copy_fixture "1_sphere/vio_A.rxt" "vio_A.rxt"
copy_fixture "1_sphere/vio_B.rxt" "vio_B.rxt"
echo "not an rxt file" > "$TEST_TMP/notes.txt"

(
    cd "$TEST_TMP"
    run_file_case basic_A.rxt basic_A.txt
    run_file_case vio_A explicit_vio_A.txt -o explicit_vio_A.txt
    "$IAD_EXECUTABLE" -M 0 -q 4 vio_A.rxt vio_B.rxt notes.txt > multiple.out 2>&1
)
assert_exists "$TEST_TMP/vio_A.txt"
assert_exists "$TEST_TMP/vio_B.txt"
assert_not_exists "$TEST_TMP/notes.txt.txt"

rm -f "$TEST_TMP/vio_A.txt" "$TEST_TMP/vio_B.txt"
(
    cd "$TEST_TMP"
    "$IAD_EXECUTABLE" -M 0 -q 4 * > glob.out 2>&1
)
assert_exists "$TEST_TMP/basic_A.txt"
assert_exists "$TEST_TMP/vio_A.txt"
assert_exists "$TEST_TMP/vio_B.txt"
assert_not_exists "$TEST_TMP/notes.txt.txt"
