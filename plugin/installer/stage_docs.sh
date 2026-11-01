#!/usr/bin/env bash
# Writes the documents every installer ships into <dest>: Readme.txt (the README's Install and Using it sections, as plain text) and
# licenses/ (the third-party licences listed in THIRD_PARTY_NOTICES.md, the notices themselves, and the project's own
# LICENSE once there is one). Called by the macOS, Windows and Linux packaging steps.
#
# Usage: installer/stage_docs.sh <dest> [version]

set -euo pipefail

DEST="${1:?Usage: stage_docs.sh <dest> [version]}"
VERSION="${2:-}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)" # plugin/
REPO="$(cd "$ROOT/.." && pwd)"                       # the README and LICENSE are at the repository root

mkdir -p "$DEST/licenses"

# Install and Using it, from their headings up to the next other one, with Markdown emphasis and links flattened.
{
    echo "REV-246${VERSION:+ $VERSION}"
    echo
    awk '/^## (Install|Using it)/{on=1; print; next} /^## /{on=0} on' "$REPO/README.md" |
        perl -0pe 's/\*\*(.+?)\*\*/$1/gs; s/\*(.+?)\*/$1/gs; s/\[([^]]+)\]\([^)]+\)/$1/g'
    echo
    echo "REV-246 is licensed as described in licenses/REV246-LICENSE.txt, if present. Source: https://github.com/tehguitarist/REV-246"
    echo
    echo "Third-party software: see the licenses folder next to this file."
} > "$DEST/Readme.txt"

cp "$ROOT/THIRD_PARTY_NOTICES.md" "$DEST/licenses/THIRD_PARTY_NOTICES.txt"
[ -f "$REPO/LICENSE" ] && cp "$REPO/LICENSE" "$DEST/licenses/REV246-LICENSE.txt"
echo "Staged docs in $DEST"
