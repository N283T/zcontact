#!/bin/sh
set -eu

if [ "$#" -ne 1 ]; then
    echo "usage: $0 STRUCTURE.pdb[.gz]|STRUCTURE.cif[.gz]" >&2
    exit 2
fi

if ! command -v gemmi >/dev/null 2>&1; then
    echo "error: gemmi is required" >&2
    exit 2
fi

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
zcontact=${ZCONTACT:-"$root/zig-out/bin/zcontact"}
if [ ! -x "$zcontact" ]; then
    echo "error: $zcontact not found; run 'zig build -Doptimize=ReleaseFast'" >&2
    exit 2
fi

tmp=${TMPDIR:-/tmp}/zcontact-gemmi-$$
trap 'rm -f "$tmp"' EXIT HUP INT TERM

"$zcontact" --mode atom --cutoff 4.0 "$1" >"$tmp"
zcontact_count=$(awk 'NR > 1 { n++ } END { print n + 0 }' "$tmp")
gemmi_count=$(gemmi contact -d 4 --ignore=1 --nosym --noh --noligand "$1" | awk 'END { print NR + 0 }')

printf 'zcontact\t%s\nGemmi\t%s\n' "$zcontact_count" "$gemmi_count"
if [ "$zcontact_count" -ne "$gemmi_count" ]; then
    echo "counts differ (altloc policy is a common intentional cause)" >&2
    exit 1
fi
