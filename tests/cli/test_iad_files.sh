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
    batch_dir=${2:-tests/rxt}
    batch_mode=${3:-default}
    case "$batch_dir" in
        tests/rxt|tests/rxt/*) ;;
        *) fail "batch directory must be under tests/rxt: $batch_dir" ;;
    esac
    case "$batch_mode" in
        default|nomc) ;;
        *) fail "unknown batch mode: $batch_mode" ;;
    esac

    manifest="$TEST_TMP/rxt_files.txt"
    (cd "$ROOT_DIR" && find "$batch_dir" -type f -name "*.rxt" | sort) > "$manifest"
    total=$(wc -l < "$manifest" | tr -d ' ')
    count=0
    while IFS= read -r path; do
        case "$path" in
            *.rxt) ;;
            *) continue ;;
        esac

        base=$(basename "$path")
        stem=${base%.rxt}
        work="$TEST_TMP/batch_$stem"
        mkdir -p "$work"
        cp "$ROOT_DIR/$path" "$work/$base"

        opts="-q 4"
        if [ "$batch_mode" = "nomc" ]; then
            opts="-M 0 $opts"
        fi
        case "$base" in
            fairway_A.rxt|fairway_E.rxt) opts="-c 0 -M 0 -q 4" ;;
            royston9_D.rxt) opts="-e 0.005 $opts" ;;
            thinh*) opts="-c 0 $opts" ;;
            tio2_vis.rxt|royston1.rxt) opts="-M 0 -q 4" ;;
            ville1.rxt) opts="-a 0 $opts" ;;
        esac

        count=$((count + 1))
        printf '[%s/%s] iad %s %s\n' "$count" "$total" "$opts" "$path" >&2
        (
            cd "$work"
            # shellcheck disable=SC2086
            "$IAD_EXECUTABLE" $opts "$base" > stdout
        ) || {
            cat "$work/stdout" >&2
            fail "iad batch case failed for $path"
        }
        assert_exists "$work/$stem.txt"
    done < "$manifest"
    echo "processed $count .rxt fixtures" >&2
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
