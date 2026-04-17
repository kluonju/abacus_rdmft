# AGENTS.md

## Cursor Cloud specific instructions

ABACUS is a C++ DFT (density functional theory) electronic structure simulation package. It builds a single binary (`abacus`) via CMake.

### Building

The build uses CMake. The standard development build command is:

```bash
cmake -B build -DBUILD_TESTING=ON -DENABLE_LIBXC=ON -DCMAKE_CXX_COMPILER=/usr/bin/g++
cmake --build build -j$(nproc)
sudo cmake --install build
```

**Gotcha:** The default `c++` alternative may point to `clang++`, which fails to link with `-lstdc++`. Always pass `-DCMAKE_CXX_COMPILER=/usr/bin/g++` to CMake.

The binary name depends on features enabled (e.g. `abacus_2p` for LCAO+MPI). A symlink `abacus` is created at install time under `/usr/local/bin/`.

### Building with Intel oneAPI

Intel oneAPI 2025.3 (icx/icpx/ifx) is installed. Source the environment first, then build with MKL:

```bash
source /opt/intel/oneapi/setvars.sh
cmake -B build_intel -DCMAKE_CXX_COMPILER=mpiicpx -DBUILD_TESTING=ON -DENABLE_LIBXC=ON
cmake --build build_intel -j$(nproc)
```

This auto-detects MKL (replacing FFTW+OpenBLAS+ScaLAPACK) and uses Intel MPI. The `setvars.sh` source is already in `~/.bashrc` for persistence.

**Gotcha:** When running Intel MPI in a container without a network fabric, set `export I_MPI_FABRICS=shm` to avoid fabric initialization errors.

### Running unit tests

```bash
cd /workspace/build                                         # or build_intel
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 OMPI_MCA_btl_vader_single_copy_mechanism=none  # OpenMPI only
export I_MPI_FABRICS=shm                                    # Intel MPI only
ctest -R "MODULE_BASE" --timeout 120 --output-on-failure    # subset
ctest --timeout 300 --output-on-failure                     # all 283 tests
```

**Gotcha:** OpenMPI refuses to run as root by default. The three `OMPI_*` env vars above are required in CI/container environments. For Intel MPI, use `I_MPI_FABRICS=shm` instead.

### Running ABACUS

To run a calculation, create a directory with `INPUT`, `STRU`, `KPT` files and pseudopotential files, then:

```bash
OMP_NUM_THREADS=2 mpirun -np 2 abacus
```

Example inputs are in `examples/`. Pseudopotentials and orbital files are in `tests/PP_ORB/`. When using examples, update `pseudo_dir` in `INPUT` to an absolute path (e.g. `/workspace/tests/PP_ORB`).

### Linting

Code formatting uses `clang-format` with the `.clang-format` config in the repo root. Check formatting with:

```bash
clang-format --dry-run --Werror <file>
```

Formatting is also enforced by pre-commit.ci after push. Local pre-commit hooks are optional per `docs/CONTRIBUTING.md`.

### Key references

- Build/test instructions: `docs/CONTRIBUTING.md`
- Integration tests: `tests/integrate/README.md`
- Examples: `examples/`
