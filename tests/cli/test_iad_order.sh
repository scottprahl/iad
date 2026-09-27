#!/bin/sh

# A row must invert the same way wherever it appears.
#
# The adaptive grid is refined around the measurement it is built for, so a
# grid built for one row is the wrong grid for the next.  It used to be built
# once per file and reused for every row after, which made the answer for a
# given wavelength depend on which wavelengths had been processed before it.
# A third of combo_1 changed depending on where the run started, and a row
# that failed in a batch could invert cleanly when rerun on its own -- so the
# obvious way to check a suspicious row silently gave a different answer than
# the batch that produced it.
#
# The test splits a file in two and requires the union of the halves to match
# the whole file exactly, row for row.

. "$(dirname "$0")/lib.sh"

announce "row order does not change the answer"

order_src="$ROOT_DIR/tests/rxt/1_sphere/homa2.rxt"
cp "$order_src" "$TEST_TMP/order.rxt"

data_rows() {
    grep -E '^[[:space:]]*[0-9]' "$1" | sort -g
}

# the fixture's data rows run from 530 to 580, so the split has to fall inside
# that range -- a split outside it puts every row in one half and the test
# passes without comparing anything
"$IAD_EXECUTABLE" -M 0 -o "$TEST_TMP/order_whole.txt" "$TEST_TMP/order.rxt" >/dev/null 2>&1 ||
    fail "iad failed on the whole file"
"$IAD_EXECUTABLE" -M 0 -l '0 552' -o "$TEST_TMP/order_lo.txt" "$TEST_TMP/order.rxt" >/dev/null 2>&1 ||
    fail "iad failed on the lower half"
"$IAD_EXECUTABLE" -M 0 -l '552.0001 100000' -o "$TEST_TMP/order_hi.txt" "$TEST_TMP/order.rxt" >/dev/null 2>&1 ||
    fail "iad failed on the upper half"

data_rows "$TEST_TMP/order_whole.txt" > "$TEST_TMP/order_whole.rows"
cat "$TEST_TMP/order_lo.txt" "$TEST_TMP/order_hi.txt" |
    grep -E '^[[:space:]]*[0-9]' | sort -g > "$TEST_TMP/order_pieces.rows"

whole_count=$(wc -l < "$TEST_TMP/order_whole.rows" | tr -d ' ')
pieces_count=$(wc -l < "$TEST_TMP/order_pieces.rows" | tr -d ' ')

lo_count=$(grep -cE '^[[:space:]]*[0-9]' "$TEST_TMP/order_lo.txt")
hi_count=$(grep -cE '^[[:space:]]*[0-9]' "$TEST_TMP/order_hi.txt")

if [ "$whole_count" -lt 10 ]; then
    fail "expected a multi-row fixture, got $whole_count rows"
fi

# a split that leaves one side empty compares nothing at all
if [ "$lo_count" -lt 4 ] || [ "$hi_count" -lt 4 ]; then
    fail "the split is lopsided: $lo_count rows below, $hi_count above"
fi

if [ "$whole_count" != "$pieces_count" ]; then
    fail "the two halves do not cover the file: whole $whole_count, pieces $pieces_count"
fi

if ! diff "$TEST_TMP/order_whole.rows" "$TEST_TMP/order_pieces.rows" > "$TEST_TMP/order.diff" 2>&1; then
    differing=$(grep -c '^<' "$TEST_TMP/order.diff")
    fail "splitting the file changed $differing of $whole_count rows
$(head -6 "$TEST_TMP/order.diff")"
fi
