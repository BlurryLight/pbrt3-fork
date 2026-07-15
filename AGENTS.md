# AGENTS.md

## Build

**Build system**: CMake + Ninja. Output directory is `build/cmake-{platform}-{buildtype}`.
When invoked by an agent, append `-agent` suffix: `build/cmake-{platform}-{buildtype}-agent`.

### Quick build (recommended)

The build script applies submodule patches then configures, builds, and tests:

```bash
bash .github/scripts/build.sh
```

### Setup (one-time)

After `git submodule update --init --recursive` or `git checkout`, the vendored
submodules need mingw-compatibility patches applied. Run once:

```bash
bash .github/scripts/setup.sh
```

Patches are idempotent and `#ifdef __MINGW32__`-guarded — safe on all platforms.
Once applied, you can invoke cmake directly (see below).

### Manual cmake

```bash
# One-time: init submodules + apply patches
git submodule update --init --recursive
bash .github/scripts/setup.sh

# Choose platform and build type
platform=$(uname -s | tr '[:upper:]' '[:lower:]')
buildtype=Release
builddir="build/cmake-${platform}-${buildtype}"

cmake -B "$builddir" -G Ninja -DCMAKE_BUILD_TYPE="$buildtype" -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build "$builddir" --parallel
```

**Agent builds** use `-agent` suffix on the build dir:
```bash
bash .github/scripts/setup.sh
builddir="build/cmake-${platform}-${buildtype}-agent"
cmake -B "$builddir" -G Ninja -DCMAKE_BUILD_TYPE="$buildtype" -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build "$builddir" --parallel
```

**CMake compatibility**: the project sets `CMAKE_MINIMUM_REQUIRED(VERSION 3.1.0)` which errors on modern CMake. Pass `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` to work around it (included above).

**Submodules are required**: cmake checks for them explicitly and aborts with a clear message if they are missing.

**Key CMake options**: `PBRT_FLOAT_AS_DOUBLE` (64-bit floats) and `PBRT_SAMPLED_SPECTRUM` (SampledSpectrum over RGBSpectrum). Both default to OFF.

**Source files are discovered via `file(GLOB ...)`** — adding new `.cpp` files requires a cmake reconfigure to pick them up.

## Test

```bash
./build/cmake-${platform}-${buildtype}/pbrt_test
```

A single binary; all tests are linked into one executable using Google Test (bundled under `src/tests/gtest/`). No per-test-file selection support out of the box — to run specific tests, use the standard Google Test `--gtest_filter` flag.

## Architecture

- `src/main/pbrt.cpp` — entrypoint for the main renderer
- `src/core/` — core types and infrastructure (geometry, spectrum, parser, API, sampling, etc.)
- `src/accelerators/` — BVH and kd-tree acceleration structures
- `src/integrators/` — rendering integrators (path, BDPT, MLT, SPPM, whitted, etc.)
- `src/materials/` — material models (Disney, glass, metal, subsurface, etc.)
- `src/lights/`, `src/cameras/`, `src/shapes/`, `src/textures/`, `src/samplers/`, `src/filters/`, `src/media/` — pluggable component directories
- `src/tools/` — utilities: `obj2pbrt`, `imgtool`, `bsdftest`, `cyhair2pbrt`
- `src/ext/` — vendored code (Hosek sky model, lodepng, rply, targa)
- `src/ext/<openexr,glog,ptex,zlib>/` — git submodules

## Code style

- C++11, `.clang-format` at `src/.clang-format`: Google style, 4-space indent
- No compiler-specific lint, format, or typecheck scripts in the repo
- Builds cleanly with Apple Clang and GCC; MSVC also supported with `/utf-8` flag

## CI

- **GitHub Actions** (`.github/workflows/cmake.yml`): builds with Ninja on `ubuntu-22.04` and `windows-latest` (msys2 ucrt64), runs `pbrt_test` on both
- **Legacy**: `.travis.yml` and `appveyor.yml` (no longer active)
