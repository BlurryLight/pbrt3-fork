#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# 1. Apply submodule patches
bash "$SCRIPT_DIR/setup.sh"

# 2. Build
WORKSPACE="$(cd "$SCRIPT_DIR/../.." && pwd)"
cd "$WORKSPACE"

build_type="${BUILD_TYPE:-Release}"
builddir="build/cmake-$build_type"

echo "=== Build type: $build_type ==="
echo "=== Build dir : $builddir ==="

launcher=()
if [ -n "${GITHUB_ACTIONS:-}" ] && command -v ccache &>/dev/null; then
    launcher=(-DCMAKE_C_COMPILER_LAUNCHER=ccache -DCMAKE_CXX_COMPILER_LAUNCHER=ccache)
fi

cmake -B "$builddir" \
    -G Ninja \
    -DCMAKE_BUILD_TYPE="$build_type" \
    -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
    "${launcher[@]}"

cmake --build "$builddir" --parallel

# 3. Test
echo "=== Running tests ==="
test_exe="$builddir/pbrt_test"
[ -f "$test_exe.exe" ] && test_exe="$test_exe.exe"
"$test_exe"
