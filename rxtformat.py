#!/usr/bin/env python3
"""Normalize the header of an iad .rxt file.

The seventeen experimental parameters at the top of an .rxt file are read by
position, not by name, so the comment beside each one is documentation and
nothing else.  Over the years those comments drifted: the same quantity is
described three different ways across the corpus, and one of the variants
says "[mm] Number of measurements" for a plain count.

This rewrites the descriptions to one wording and lines the '#' up in column
12, without touching the numbers.

What is preserved:

* the comment block between the IAD1 line and the first parameter
* standalone comment lines inside the header, kept above the parameter they
  preceded
* any text one of the first six parameters carried beyond its standard
  description -- the literature citations for refractive indices, the note
  about a petri dish

What is not preserved: the comments on the sphere count and the ten sphere
parameters.  Those are replaced outright with the standard wording, since
old files swap the port names between blocks and nothing in them is worth
keeping.
* blank lines, as the grouping used by the reference files
* the whole data section, verbatim, including the column-letter line
* the line ending the file already used

Numbers are parsed to decide what the file describes -- whether spheres are
in use, whether the exit port is closed -- but are written back exactly as
they appeared, so no value is ever re-rounded.

Usage:
    rxtformat.py FILE...             write the result to stdout
    rxtformat.py --in-place FILE...  rewrite each file
    rxtformat.py --check FILE...     report which files would change
"""

import argparse
import os
import re
import sys

# Each header slot: the description to write, and the patterns that count as
# an existing description of that slot.  For the first six slots anything in
# the old comment that is not matched by one of the patterns is kept and
# appended; the sphere count and the two sphere blocks are always rewritten
# from scratch (see FIRST_SPHERE_SLOT).
NUMBER = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"

SAMPLE_INDEX = "Index of refraction of the sample"
SLIDE_INDEX = "Index of refraction of the top and bottom slides"
SAMPLE_THICKNESS = "[mm] Thickness of sample"
SLIDE_THICKNESS = "[mm] Thickness of slides"
BEAM = "[mm] Diameter of illumination beam"
RSTD = "Reflectance of the calibration standard"
NUM_SPHERES = "Number of spheres used at a time to make a measurement"
SPHERE_D = "[mm] Sphere Diameter"
SAMPLE_PORT = "[mm] Sample Port Diameter"
ENTRANCE_PORT = "[mm] Light Entrance Port Diameter"
EXIT_PORT = "[mm] Light Exit Port Diameter"
DETECTOR_PORT = "[mm] Detector Port Diameter"
WALL = "Reflectance of the sphere wall"

# Header slot of the sphere count; it and every slot after it describe the
# spheres.
FIRST_SPHERE_SLOT = 6

REFLECTION_BLOCK = "Sphere configuration for reflectance measurements"
TRANSMISSION_BLOCK = "Sphere configuration for transmittance measurements"

# Patterns are matched case-insensitively against the old comment with the
# leading '#' and whitespace removed.  Order matters only in that the first
# match wins, so the more specific spelling comes first.
KNOWN = {
    SAMPLE_INDEX: [r"index of refraction of (?:the )?sample"],
    SLIDE_INDEX: [
        r"index of refraction of (?:the )?top and bottom slides?",
        r"index of refraction of (?:the )?(?:top|bottom|both) slides?",
        r"index of refraction of (?:the )?slides?",
    ],
    SAMPLE_THICKNESS: [r"(?:\[mm\]\s*)?thickness of (?:the )?sample"],
    SLIDE_THICKNESS: [r"(?:\[mm\]\s*)?thickness of (?:the )?slides?"],
    BEAM: [
        r"(?:\[mm\]\s*)?diameter of (?:the )?illumination beam",
        r"(?:\[mm\]\s*)?(?:illumination )?beam diameter",
    ],
    RSTD: [
        r"reflect\w* of (?:the )?(?:reflectance )?(?:calibration )?standard",
    ],
}

# Read_Data_Line in iad_io.w takes the fixed columns (M_R, M_T, M_U, r_w,
# r_w_t, rstd_r, rstd_t) in this order and stops at the count given, so the
# count alone says what the columns are.  A leading wavelength column is
# optional and is recognised at read time by its first value being greater
# than one, which is why it is not part of this list.
#
# The label each fixed column would carry in the newer format, taken from
# Read_Data_Line_Per_Labels.  The correspondence is exact: 'w' sets both wall
# reflectances just as the fourth fixed column does, 'W' sets only the
# transmission one, and 'R' and 'T' are the two calibration standards.  A
# count is therefore rewritten as labels without changing what iad reads --
# and the file then says what its columns are instead of leaving the reader
# to remember an order.
FIXED_LETTERS = ["r", "t", "u", "w", "W", "R", "T"]
LAMBDA_LETTER = "L"

# What every label means, from Read_Data_Line_Per_Labels, named as in the
# column-label table of docs/manual.tex.  A label line is written with this
# spelled out above it, so a reader does not have to know the table by heart.
LETTER_NAMES = {
    "a": ("albedo", ""), "A": ("mu_a", ""),
    "b": ("opt depth", ""), "B": ("beam diam", ""),
    "c": ("unscat R", ""), "C": ("unscat T", ""),
    "d": ("d_sample", ""), "D": ("d_slide", ""),
    "e": ("tolerance", ""), "E": ("b_slide", ""),
    "f": ("wall first", ""),
    "F": ("mu_s", ""), "g": ("anisotropy", ""),
    "i": ("incid angle", ""), "L": ("lambda", ""),
    "M": ("spheres", ""), "n": ("n_sample", ""),
    "N": ("n_slide", ""), "q": ("quad pts", ""),
    "r": ("M_R", ""), "R": ("rstd", ""),
    "S": ("spheres", ""), "t": ("M_T", ""),
    "T": ("tstd", ""), "u": ("M_U", ""),
    "w": ("r_w", ""), "W": ("t_w", ""),
}


def names_for_letters(letters):
    """Two-part name for each label, unrecognised letters standing for
    themselves.  Every part is short enough for one field, so nothing has to
    be chopped mid-word."""
    return [LETTER_NAMES.get(letter, (letter, "")) for letter in letters]


def spell_out(pairs):
    """One-line rendering of the column names, for messages and comparisons."""
    return ", ".join(" ".join(part for part in pair if part) for pair in pairs)

# Every letter Read_Data_Line_Per_Labels understands.  Anything else on a
# label line is left alone rather than guessed at.
COLUMN_LETTERS = set("aAbBcCdDeEfFgiLMnNqrRStTuwW")

LABEL_WIDTH = 10
FIELD_WIDTH = 11

BLOCK_PATTERNS = [
    r"sphere configuration for \w+ measurements(?:\s*\(unused\))?",
    r"properties of (?:the )?(?:\w+ )?sphere[^#]*",
]

COMMENT_COLUMN = 11
FIRST_LINE = "IAD1".ljust(COMMENT_COLUMN) + "# required first four characters"


# The old IAD format (first line "IAD" rather than "IAD1") has three more
# header values after the sphere blocks: the fraction of light that hits the
# sphere wall first, whether M_R includes the specular reflection, and whether
# M_T includes the unscattered transmission.  IAD1 has no place for them in
# the header, so they become three constant data columns with the labels that
# set the same quantities, and nothing the old file said is lost.
LEGACY_LETTERS = ["f", "c", "C"]


class RxtError(Exception):
    """Raised when a file cannot be understood as an .rxt header."""


def split_comment(line):
    """Return (code, comment) for one line; comment excludes the '#'."""
    hash_at = line.find("#")
    if hash_at < 0:
        return line, None
    return line[:hash_at], line[hash_at + 1:].strip()


def strip_known(comment, description):
    """Remove a recognised description, returning whatever text is left.

    Args:
        comment: The old comment text, without the leading '#'.
        description: The standard description this slot should carry.

    Returns:
        The remainder of the comment once a recognised spelling of the
        description has been removed.  A comment that matches nothing is
        returned whole, so unfamiliar notes are never dropped.
    """
    if comment is None:
        return ""

    text = comment.strip()
    for pattern in KNOWN.get(description, []):
        match = re.match(pattern, text, flags=re.IGNORECASE)
        if match:
            return text[match.end():].strip()
    return text


def describe(description, extra):
    """Join a standard description with any preserved text."""
    if not extra:
        return description
    if extra.startswith(",") or extra.startswith(";"):
        return description + extra
    return description + " " + extra


def parameter_line(value_text, description, extra=""):
    """Format one parameter with the '#' in column 12."""
    text = describe(description, extra)
    return ("%-*s# %s" % (COMMENT_COLUMN, value_text, text)).rstrip()


def comment_line(text):
    """Format a standalone comment with the '#' in column 12."""
    return ("%-*s# %s" % (COMMENT_COLUMN, "", text)).rstrip()


COLUMN_WORDS = ("m_r", "m_t", "m_u", "r_w", "t_w", "rstd", "tstd", "lambda",
                "nm", "wavelength", "mr", "mt", "mu")


def looks_like_column_list(text):
    """True when text is just a list of column names, not a real note."""
    if not text:
        return False
    words = [w for w in re.split(r"[\s,;]+", text.strip().lower()) if w]
    if not words:
        return False
    return all(any(w.startswith(c) for c in COLUMN_WORDS) for w in words)


def drop_old_header(lines, data_line):
    """Skip a column-header comment sitting just above the data."""
    index = data_line
    while index < len(lines):
        code, comment = split_comment(lines[index])
        if code.strip():
            break
        if comment is None:
            break
        if looks_like_column_list(comment) or looks_like_old_parameter(comment):
            index += 1
            continue
        break
    return index


def looks_like_old_parameter(text):
    """True for a commented-out header line such as '0.000  # Fraction ...'.

    Old files kept parameters iad no longer reads by commenting them out;
    they describe nothing about the file as it now stands.
    """
    return re.match(NUMBER + r"\s*(?:#|$)", text.strip()) is not None


def first_data_value(lines, data_line, shared):
    """The first number of the first data row, or None if there is none.

    Read_Data_Line decides whether a row begins with a wavelength by looking
    at this one value: greater than one means a wavelength, anything else is
    already M_R.  The column count cannot be used instead, because a file may
    carry extra columns for other reasons.
    """
    if shared:
        return float(shared[0]["text"])
    for line in lines[data_line:]:
        code, _ = split_comment(line)
        pieces = code.split()
        if pieces:
            try:
                return float(pieces[0])
            except ValueError:
                return None
    return None


def count_data_columns(lines, data_line):
    """Number of whitespace-separated values on the first real data row."""
    for line in lines[data_line:]:
        code, _ = split_comment(line)
        pieces = code.split()
        if pieces:
            return len(pieces)
    return 0


def header_comment(names):
    """A commented column-header line lined up with the data beneath it."""
    cells = "".join("%-*s" % (LABEL_WIDTH, name) for name in names)
    return ("#" + cells).rstrip()


def label_line(letters):
    """Write the column labels one per fixed-width, tab-separated field."""
    return "\t".join("%-*s" % (FIELD_WIDTH, letter) for letter in letters)


def fit_number(text):
    """Shorten a number to the field width, rounding only if it must.

    A value already short enough is passed through untouched, so nothing is
    re-rounded for the sake of it.  Longer values lose digits: this is the one
    place the program changes the data rather than the documentation.
    """
    if len(text) <= FIELD_WIDTH:
        return text
    try:
        value = float(text)
    except ValueError:
        return text
    # fixed point keeps the trailing zeros that %g strips, so 0.0173065998
    # stays 0.017306600 instead of looking less precise than 1.730660e-2
    whole = len("%d" % abs(value)) + (1 if value < 0 else 0)
    decimals = FIELD_WIDTH - whole - 1
    best = "%.*f" % (decimals, value) if decimals > 0 else None
    for digits in range(FIELD_WIDTH, 0, -1):
        for candidate in ("%.*g" % (digits, value), compact_exponent(value, digits)):
            if len(candidate) <= FIELD_WIDTH:
                if best is None or significant_digits(candidate) > significant_digits(best):
                    best = candidate
        if best is not None and significant_digits(best) >= digits:
            break
    return best if best is not None else "%.0e" % value


def compact_exponent(value, digits):
    """Exponential form with the exponent trimmed, as in 1.4e-4.

    A small value written plainly loses almost everything at this width:
    -0.00014358 becomes -0.0001, one significant figure.  Dropping the
    padding zeros from the exponent buys back a digit or two.
    """
    text = "%.*e" % (digits - 1, value)
    mantissa, _, exponent = text.partition("e")
    sign = "-" if exponent[0] == "-" else ""
    return mantissa + "e" + sign + exponent[1:].lstrip("0").rjust(1, "0")


def significant_digits(text):
    """Count the significant figures a rendering actually carries."""
    body = text.split("e")[0].lstrip("-+")
    return len(body.replace(".", "").lstrip("0")) or 1


def data_row(values, comment):
    """One row of data, each value in its own field, comment kept at the end."""
    cells = "\t".join("%-*s" % (FIELD_WIDTH, fit_number(v)) for v in values)
    if comment:
        return cells + "\t# " + comment
    return cells


def name_header(names):
    """Two comment lines naming the columns, aligned with the data below.

    Long names such as "sample thickness" are unreadable squeezed into seven
    characters, so each is allowed a second line underneath.  The '#' eats one
    character of the first field so the columns still line up.
    """
    tops = [pair[0] for pair in names]
    bottoms = [pair[1] for pair in names]

    def render(cells):
        first = "#%-*s" % (FIELD_WIDTH - 1, cells[0][:FIELD_WIDTH - 1])
        rest = "\t".join("%-*s" % (FIELD_WIDTH, c) for c in cells[1:])
        return (first + "\t" + rest).rstrip() if rest else first.rstrip()

    lines_out = [render(tops)]
    if any(bottoms):
        lines_out.append(render(bottoms))
    return lines_out


def tokenize(lines):
    """Collect the numbers and column letters in file order.

    Returns:
        A list of dicts with the token text, the line it came from, that
        line's comment, and whether the token is the last one on its line.
    """
    tokens = []
    for index, line in enumerate(lines):
        code, comment = split_comment(line)
        if index == 0:
            code = re.sub(r"IAD1?", " ", code, count=1)
        pieces = code.replace(",", " ").split()
        for position, piece in enumerate(pieces):
            tokens.append({
                "text": piece,
                "line": index,
                "comment": comment,
                "last_on_line": position == len(pieces) - 1,
            })
    return tokens


def is_number(text):
    """True when the token is a number rather than a column letter."""
    return re.fullmatch(NUMBER, text) is not None


def leading_comments(lines, first_param_line):
    """Comment lines between the IAD1 line and the first parameter."""
    kept = []
    for line in lines[1:first_param_line]:
        code, comment = split_comment(line)
        if code.strip():
            break
        if comment is not None:
            kept.append(comment)
    return kept


def interleaved_comments(lines, start, stop):
    """Standalone comment lines strictly between two parameter lines."""
    kept = []
    for line in lines[start + 1:stop]:
        code, comment = split_comment(line)
        if code.strip():
            continue
        if comment is not None and not is_block_header(comment):
            kept.append(comment)
    return kept


def is_block_header(comment):
    """True for the 'Sphere configuration ...' headings we regenerate."""
    text = comment.strip()
    for pattern in BLOCK_PATTERNS:
        if re.match(pattern, text, flags=re.IGNORECASE):
            return True
    return False


def format_data(lines_in, extra=()):
    """Lay the data out in fixed-width, tab-separated fields.

    Blank lines and whole-line comments are passed through; a row of numbers
    is rewritten, with the values in extra appended, and any comment trailing
    it is kept.
    """
    out = []
    for line in lines_in:
        code, comment = split_comment(line)
        values = code.split()
        if not values:
            out.append(line.rstrip())
            continue
        out.append(data_row(values + list(extra), comment))
    return out


def format_rxt(text, warn=lambda message: None):
    """Return the file text with a normalized header.

    Args:
        text: Entire contents of an .rxt file.
        warn: Called with a message for anything accepted but not quite right.

    Returns:
        The reformatted contents.

    Raises:
        RxtError: The file does not begin with IAD, has too few numbers, or
            has a measurement count that is not 1 to 7.
    """
    newline = "\r\n" if "\r\n" in text else "\n"
    lines = text.replace("\r\n", "\n").split("\n")

    if not lines or "IAD" not in lines[0]:
        raise RxtError("does not start with IAD1")
    old_format = "IAD1" not in lines[0]

    tokens = tokenize(lines)
    # positions are carried along rather than looked up later: two tokens on
    # the same line with the same text are equal as dicts, so searching for
    # one finds the earlier one.  A file written entirely on one line, like
    # terse_A.rxt, is all duplicates.
    numbered = [(i, t) for i, t in enumerate(tokens) if is_number(t["text"])]
    if len(numbered) < 17:
        raise RxtError("expected 17 header numbers, found %d" % len(numbered))

    header = [t for _, t in numbered[:17]]
    header_count = 17
    legacy_values = []
    if old_format:
        header_count = 20
        if len(numbered) < 20:
            raise RxtError("old IAD format needs 20 header numbers, found %d"
                           % len(numbered))
        legacy_values = [t["text"] for _, t in numbered[17:20]]
        warn("old IAD format converted to IAD1; its three extra header "
             "values are now the constant columns f, c, and C")
    last_header = numbered[header_count - 1][1]
    after = tokens[numbered[header_count - 1][0] + 1:]

    # The count of measurements is a number; a column-letter line is not.
    count_token = None
    if after and is_number(after[0]["text"]):
        count_token = after[0]
        count = float(count_token["text"])
        if count != int(count) or not 1 <= count <= 7:
            raise RxtError("expected the number of measurements (1-7) after "
                           "the header, found %s" % count_token["text"])

    num_spheres = float(header[6]["text"])
    unused = " (unused)" if num_spheres == 0 else ""

    # the IAD1 line is followed by exactly one blank line, then the comment
    # block that names the file, then one blank line before the parameters
    out = [FIRST_LINE, ""]
    # the block that names the file stays flush left, the way it is written
    # by hand; only the per-parameter comments get lined up
    preamble = leading_comments(lines, header[0]["line"])
    for comment in preamble:
        out.append(("# " + comment).rstrip())
    if preamble:
        out.append("")

    slots = [
        SAMPLE_INDEX, SLIDE_INDEX, SAMPLE_THICKNESS, SLIDE_THICKNESS,
        BEAM, RSTD, NUM_SPHERES,
        SPHERE_D, SAMPLE_PORT, ENTRANCE_PORT, DETECTOR_PORT, WALL,
        SPHERE_D, SAMPLE_PORT, EXIT_PORT, DETECTOR_PORT, WALL,
    ]

    for slot in range(17):
        token = header[slot]
        description = slots[slot]

        # blank line and heading before each sphere block, and before the
        # sphere count, matching the layout of the reference files
        if slot == 6:
            out.append("")
        if slot == 7:
            out.append("")
            out.append(comment_line(REFLECTION_BLOCK + unused))
        if slot == 12:
            out.append("")
            out.append(comment_line(TRANSMISSION_BLOCK + unused))

        # standalone comments inside the sphere section are old block
        # headings ("Reflection Sphere" and the like) that the headings
        # above replace, so only those among the first six are kept
        if 0 < slot < FIRST_SPHERE_SLOT:
            for comment in interleaved_comments(lines, header[slot - 1]["line"],
                                                token["line"]):
                out.append(comment_line(comment))

        # the sphere comments are the ones that drifted most, with port
        # names swapped between blocks and earlier runs of this script
        # appended, so nothing of the old text is worth keeping
        if slot >= FIRST_SPHERE_SLOT:
            extra = ""
        else:
            comment = token["comment"] if token["last_on_line"] else None
            extra = strip_known(comment, description)

        if description is EXIT_PORT and float(token["text"]) == 0:
            extra = "(closed)"

        if description is SLIDE_INDEX and float(token["text"]) == 1:
            if "(none)" not in extra:
                extra = ("(none) " + extra).strip()

        value = token["text"]
        if description is SLIDE_THICKNESS:
            # a slide with the index of air is no slide at all, so whatever
            # thickness was given for it describes nothing
            if float(header[1]["text"]) == 1 and float(value) != 0:
                warn("slide index is 1.0, so slide thickness %s set to 0" % value)
                value = "0"
            if float(value) == 0 and "(none)" not in extra:
                extra = ("(none) " + extra).strip()

        out.append(parameter_line(value, description, extra))

    if count_token is not None:
        data_tokens = after[1:]
        header_ends_on = count_token["line"]
    else:
        data_tokens = after
        header_ends_on = last_header["line"]

    # Usually the data starts on its own line and can be copied across
    # untouched.  It does not have to: a file may put the header and the data
    # on one line, and then the leftover tokens are all that is left of it.
    shared = [t for t in data_tokens if t["line"] == header_ends_on]
    if shared:
        data_line = header_ends_on + 1
    else:
        data_line = data_tokens[0]["line"] if data_tokens else len(lines)

    wrote_labels = False

    if count_token is not None:
        # the count says which fixed columns are present; a leading wavelength
        # column is optional, so it is inferred by counting the data itself
        count = int(float(count_token["text"]))
        leading = first_data_value(lines, data_line, shared)
        has_lambda = leading is not None and leading > 1

        # the comment on the count is always some description of the
        # columns ("Two measurements, i.e., M_R & M_T"), often out of date;
        # the label line written below replaces it
        letters = FIXED_LETTERS[:min(count, len(FIXED_LETTERS))]
        if has_lambda:
            letters = [LAMBDA_LETTER] + letters
        if legacy_values:
            letters = letters + LEGACY_LETTERS

        out.append("")
        out.extend(name_header(names_for_letters(letters)))
        out.append(label_line(letters))
        wrote_labels = True

        # a column header already in the file is replaced by the one just
        # written; any other comment down there is somebody's note and stays
        data_line = drop_old_header(lines, data_line)
    else:
        # a label line names its own columns; keep the letters, fix the spacing
        label_row = data_tokens[0]["line"] if data_tokens else -1
        letters = [t["text"] for t in data_tokens if t["line"] == label_row]
        if letters and all(len(x) == 1 and x in COLUMN_LETTERS for x in letters):
            spelled = spell_out(names_for_letters(letters))
            generated_header = set(name_header(names_for_letters(letters)))
            out.append("")
            # notes written above a label line are kept; the spelled-out
            # column list is dropped because it is rewritten just below, and
            # matching it exactly is safer than guessing from its words --
            # "sample index" and "anisotropy" are column names too
            for note in lines[last_header["line"] + 1:label_row]:
                code, note_comment = split_comment(note)
                if code.strip() or note_comment is None:
                    continue
                if note_comment.strip() == spelled:
                    continue
                if ("#" + note_comment) in generated_header:
                    continue
                if any(note.rstrip() == g for g in generated_header):
                    continue
                if looks_like_column_list(note_comment):
                    continue
                if looks_like_old_parameter(note_comment):
                    continue
                out.append(comment_line(note_comment))
            if legacy_values:
                letters = letters + LEGACY_LETTERS
            out.extend(name_header(names_for_letters(letters)))
            out.append(label_line(letters))
            data_line = label_row + 1
            wrote_labels = True

    if shared:
        if not wrote_labels:
            out.append("")
        out.append(data_row([t["text"] for t in shared] + legacy_values, None))

    # the data section keeps exactly one blank line ahead of it, however many
    # the original had
    tail = format_data(lines[data_line:], legacy_values)
    while tail and not tail[0].strip():
        tail = tail[1:]
    # a column header belongs directly above its data; anything else gets a
    # blank line to separate it
    if not wrote_labels:
        out.append("")
    out.extend(tail)

    return newline.join(out)


def main():
    """Reformat the headers of the .rxt files named on the command line."""
    parser = argparse.ArgumentParser(
        description="Normalize the header comments of iad .rxt files.")
    parser.add_argument("files", nargs="+", metavar="FILE")
    parser.add_argument("--in-place", action="store_true",
                        help="rewrite each file instead of printing it")
    parser.add_argument("--check", action="store_true",
                        help="report which files would change, change nothing")
    args = parser.parse_args()

    if args.in_place and args.check:
        parser.error("--in-place and --check cannot be combined")

    would_change = 0
    failed = 0

    for path in args.files:
        if not path.endswith(".rxt"):
            print("skipping %s: not a .rxt file" % path, file=sys.stderr)
            continue
        try:
            with open(path, encoding="utf-8", newline="") as handle:
                original = handle.read()
            formatted = format_rxt(
                original,
                lambda message: print("%s: warning: %s" % (path, message),
                                      file=sys.stderr))
        except (OSError, RxtError, ValueError) as exc:
            print("%s: %s" % (path, exc), file=sys.stderr)
            failed += 1
            continue

        if args.check:
            if formatted != original:
                print("would reformat %s" % path)
                would_change += 1
        elif args.in_place:
            if formatted != original:
                with open(path, "w", encoding="utf-8", newline="") as handle:
                    handle.write(formatted)
                print("reformatted %s" % os.path.basename(path))
        else:
            sys.stdout.write(formatted)

    if failed:
        return 1
    if args.check and would_change:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
