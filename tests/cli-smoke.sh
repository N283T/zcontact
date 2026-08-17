#!/bin/sh
set -eu

bin=$1
fixture=$2
tmp=$(mktemp -d "${TMPDIR:-/tmp}/zcontact-cli.XXXXXX")
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

"$bin" --version >"$tmp/version.out" 2>"$tmp/version.err"
test ! -s "$tmp/version.err"
awk '$1 == "zcontact" && $2 == "0.1.0" { ok=1 } END { exit !ok }' "$tmp/version.out"

"$bin" "$fixture" >"$tmp/default.tsv"
awk -F '\t' 'NR == 1 && $1 == "residue_index1" { ok=1 } END { exit !ok }' "$tmp/default.tsv"

# Explicit asymmetric sides must control output orientation.
"$bin" --mode atom --scope all --select1 chain:A --select2 chain:B "$fixture" >"$tmp/oriented.tsv"
awk -F '\t' 'NR > 1 && ($3 != "A" || $12 != "B") { exit 1 } END { exit NR < 2 }' "$tmp/oriented.tsv"

if "$bin" --select '' "$fixture" >/dev/null 2>/dev/null; then
    echo "empty selection unexpectedly succeeded" >&2
    exit 1
fi

# Failures before publication must preserve an existing destination.
printf 'KEEP\n' >"$tmp/existing.tsv"
printf 'not a structure\n' >"$tmp/bad.pdb"
if "$bin" "$tmp/bad.pdb" -o "$tmp/existing.tsv" >/dev/null 2>/dev/null; then
    echo "malformed input unexpectedly succeeded" >&2
    exit 1
fi
test "$(cat "$tmp/existing.tsv")" = KEEP

# Batch --select1 alone inherits select2 exactly as the single-file CLI does.
mkdir "$tmp/in"
cp "$fixture" "$tmp/in/x.pdb"
"$bin" batch "$tmp/in" -o "$tmp/out" --select1 chain:A --quiet
awk -F '\t' '$1 == "select1" && $2 == "chain:A" { a=1 } $1 == "select2" && $2 == "chain:A" { b=1 } END { exit !(a && b) }' "$tmp/out/run.tsv"
