# hummingbird
Solving the self-adjoint angular flux transport equation using spectral elements on Cartesian geometry.

## Installation
Hummingbird uses [pixi](https://pixi.sh/) to manage its build tools and dependencies, so no manual installation of `cmake`, compilers, or third-party libraries is required.

[pixi](https://pixi.sh/) is a cross-platform package and workspace manager built on the conda ecosystem. It resolves all of the project's dependencies (compilers, `cmake`, `ninja`, and libraries like `armadillo` and `vtk`) from `conda-forge` into a self-contained environment under `.pixi/`, pinned by a lockfile ([`pixi.lock`](./pixi.lock)) so every contributor and CI run builds against identical dependency versions. This means:
- No system-wide installation of build tools or libraries, and no risk of conflicting with other projects' dependencies on your machine.
- Reproducible builds, since the lockfile guarantees everyone resolves the same dependency versions instead of "works on my machine" drift.
- A single command (`pixi install`) sets up everything needed to build, test, and generate documentation for the project.

### Installing pixi
If you don't already have pixi installed, install it with:

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

See the [pixi documentation](https://pixi.sh/latest/installation/) for other installation methods (Homebrew, Winget, etc.).

### Setting up the environment
From the repository root, run:

```bash
pixi install
```

This reads [`pixi.toml`](./pixi.toml) / [`pixi.lock`](./pixi.lock) and creates a local environment (in `.pixi/`) containing all the dependencies listed below, resolved from the `conda-forge` channel. Supported platforms are `linux-64` and `osx-arm64`.

## Usage
Hummingbird defines its build/test/documentation steps as [pixi tasks](https://pixi.sh/latest/workspace/advanced_tasks/), which run inside the pixi-managed environment.

### `pixi run dev`
Runs the full development workflow (`cmake --workflow --preset dev`), which:
1. **Configures** the project with CMake (`debug` preset).
2. **Builds** all targets, including the `test_hummingbird` test binary.
3. **Tests** the build by running the full test suite with CTest (`test-all` preset), stopping on the first failure and printing output for any failures.

```bash
pixi run dev
```

### `pixi run docs`
Builds the [Doxygen](https://www.doxygen.nl/) HTML documentation into `docs/documentation.html`:

```bash
pixi run docs
```

### Testing
Tests are written with [GoogleTest](https://github.com/google/googletest) and discovered automatically via CTest (`gtest_discover_tests`). `pixi run dev` builds and runs the entire suite in one step. If you've already configured/built the project and just want to re-run the tests:

```bash
ctest --test-dir build --output-on-failure
```

### Manual CMake configuration
For editor/IDE integration (e.g. CLion, VS Code with the CMake Tools extension), a `pixi` CMake preset is also provided, which points CMake at the compilers and tools installed by pixi (`.pixi/envs/default/bin`) instead of using `pixi run` directly:

```bash
cmake --preset pixi
```

## Dependencies
All dependencies below are installed automatically by `pixi install` (see [`pixi.toml`](./pixi.toml) for exact version constraints):

| Package | Purpose |
| --- | --- |
| [cmake](https://cmake.org/) | Build configuration |
| [ninja](https://ninja-build.org/) | Build system used by CMake |
| [armadillo](https://arma.sourceforge.net/) | Linear algebra |
| [doxygen](https://www.doxygen.nl/) | Documentation generation |
| [graphviz](https://graphviz.org/) | Diagrams for generated documentation |
| [gtest](https://github.com/google/googletest) | Unit testing framework |
| [nlohmann_json](https://github.com/nlohmann/json) | JSON parsing for simulation input files |
| [gmsh](https://gmsh.info/) | Mesh generation/import |
| [vtk](https://vtk.org/) | Visualization output |
| [fmt](https://github.com/fmtlib/fmt) | String formatting |

## Acknowledgements
This project contains code either directly copied or repurposed from [starling](https://github.com/lipohlmann/starling/), licensed under BSD 3-Clause.

### Third Party
This project contains code from [TinyExpr++](https://github.com/Blake-Madden/tinyexpr-plusplus), the C++ version of the [TinyExpr](https://github.com/codeplea/tinyexpr) library. These include [tinyexpr.h](./include/utils/tinyexpr.h) and [tinyexpr.cpp](./src/utils/tinyexpr.cpp). [TinyExpr++](https://github.com/Blake-Madden/tinyexpr-plusplus) is licensed under zlib, the language for which can be found [here](https://github.com/Blake-Madden/tinyexpr-plusplus/blob/tinyexpr%2B%2B/LICENSE). No part of this code has been modified from its original source aside from the include path in [tinyexpr.cpp](./src/utils/tinyexpr.cpp) to the header file and a reformatting by clang-format to match the project.
