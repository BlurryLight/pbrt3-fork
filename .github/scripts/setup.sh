#!/usr/bin/env bash
# setup.sh -- apply mingw patches to submodules (idempotent).
# Run once after checkout/submodule update, or let build.sh handle it.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
WORKSPACE="$(cd "$SCRIPT_DIR/../.." && pwd)"
PATCHES_DIR="$WORKSPACE/.github/patches"

apply_one() {
    local repo="$1"
    local patch="$2"

    if [ ! -f "$patch" ]; then
        echo "   [SKIP] patch file not found: $patch"
        return 0
    fi

    echo "  -> $(basename "$repo")"

    if (cd "$repo" && git apply --reverse --check "$patch" 2>/dev/null); then
        echo "     (already applied, skip)"
        return 0
    fi

    if (cd "$repo" && git apply --check "$patch" 2>/dev/null); then
        (cd "$repo" && git apply "$patch")
        echo "     applied"
        return 0
    fi

    echo "  ERROR: patch does not apply cleanly: $patch" >&2
    return 1
}

echo "=== Applying submodule patches ==="
apply_one "$WORKSPACE/src/ext/openexr" "$PATCHES_DIR/openexr-mingw.patch"
apply_one "$WORKSPACE/src/ext/ptex"     "$PATCHES_DIR/ptex-mingw.patch"
apply_one "$WORKSPACE/src/ext/glog"     "$PATCHES_DIR/glog-mingw.patch"
echo "=== Done ==="
