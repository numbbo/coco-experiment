Working with the `cocoex` code base
===================================

## Before you start

To work with the `coco-experiment` source code, you need quite a few tools in very specific versions.
We provide a Conda environment file `env.yaml` and recommend you use it to setup your development environment.
It is regularly tested on Linux, Windows and macOS and is what is used in CI.

Following these steps should setup a fresh `cocoex` environment:

1. Create or update a new conda environment with all the requirements by running
 
   ```sh
   conda env create -f env.yaml
   ``` 
   or if you've already setup a `cocoex` environment
   ```sh
   conda env update -f env.yaml --prune
   ``` 

1. Activate environment
   ```sh
   conda activate cocoex
   ```

You now have all required dependencies to work on the code base.
If you notice that something is missing, please let us know by opening an [issue](https://github.com/numbbo/coco-experiment/issues/new/choose).

### Reusing a conda environment

Say you already have a conda environment named `foo` that you want to reuse because it contains some additional packages you need for development.
Then, once you've activated the environment, you can run

```sh
conda env update --file env.yaml
```

and it will install all the required development dependencies into your existing environment.

## Working with the Experimental Code

Before any of the bindings can be built and after every change to the core files (in `src/`), you need to run `scripts/fabricate`.
The job of `fabricate` is to bundle the C++ files and place them in all the right places so that the language specific build tools can find them when building the respective bindings.
`fabricate` also updates the binding metadata to include version information and any other changes that need to be made from one build to the next.
It is harmless to run `fabricate` too often, but forgetting to run it means you might be using stale sources.

### Coding conventions for the core

The core in `src/` is C++17 that is amalgamated into a single `coco.cpp` and exposes a C interface (`coco.h`).

- **No exceptions.** COCO is compiled with `-fno-exceptions` where the compiler supports it. Report errors with
  `coco_error()`, `coco_warning()` etc. and avoid library functions that signal errors by throwing
  (such as `std::vector::at()` or `std::stoi()`).
- Manage memory with standard containers and RAII instead of `coco_allocate_memory()`/`coco_free_memory()`,
  which remain part of the C interface.
- Problems, transformations and loggers are subclasses of `coco_problem_s`. Suites are subclasses of
  `coco_suite_s` that implement `get_problem()` (and `get_instances_by_year()` if needed). Observers are
  subclasses of `coco_observer_s` that implement `add_logger()` (and `signal_restart()` if supported).
- Problems are owned through `coco_problem_ptr_t` (a `std::unique_ptr<coco_problem_t>`). Functions that wrap
  a problem take it by value and the caller hands it over with `std::move()`, e.g.
  `problem = transform_obj_shift(std::move(problem), fopt);`. Ownership is only released to a raw pointer
  at the C interface.

### Running unit tests

Change to the unit tests directory
```sh
cd test/unit-test
```

Rerun `fabricate` if not already done
```sh
python ../../scripts/fabricate
```

Build tests using [`cmake`](https://cmake.org/cmake/help/latest/manual/cmake.1.html)
```sh
cmake -B build
cmake --build build
```

Run tests using [`ctest`](https://cmake.org/cmake/help/latest/manual/ctest.1.html)
```sh
ctest --test-dir build     
```
or on Windows
```sh
ctest --test-dir build -C Debug     
```

### Golden master tests

The golden master tests check that the results of COCO (function values, constraint values and logger output) stay exactly the same.
Run them after every change that should not change results.

```sh
cd test/golden-test
python ../../scripts/fabricate
cmake -B build
cmake --build build
ctest --test-dir build
```

References are kept per platform.
See `test/golden-test/README.md` for how to create them for a new platform.

### Regression tests

Note that I'm not sure the regression tests are really useful at the moment. But here goes:

Change to the regression tests directory
```sh
cd test/regression-test
```

The regression tests require the `cocoex` Python package. Install it first
```sh
python -m pip uninstall -y coco-experiment
python -m pip install ../../build/python/  
```

Now run the regression tests

```sh
python test_suites.py   
```
and 
```sh
python test_loggers.py
```

### Integration tests

Still need to fix those up.

### `cocoex` Python package

Change to the `cocoex` Python package directory
```sh
cd build/python
```

Install the package from source
```sh
pip uninstall -y coco-experiment
pip install .
```

#### Compilation and installation hints

- Under macOS with an ARM processor it may be necessary to run `arch -arm64 pip install .` instead.

- When installing from source on older machines, it might be necessary to explicitly enable C99 mode for the C compiler.
  You can do this running `CFLAGS="-std=c99" pip install .`. 

#### Testing the `cocoex` Python package

Install and run `pytest`
```sh
pip install pytest
python -m pytest test
```

## Release a new version

Once you've decided a commit is stable and ready to be released, tag it with the appropriate version and push the tag to Github:

```sh
git tag -a -m 'Version X.Y.Z' vX.Y.Z
git push --tags
```

Note: You can use suffixes for pre-releases such as version 2.6.99rc4 (release candidate 4) or 2.6.99b1 (beta 1).

Pushing a tag starting with `v` triggers the `tag-release` workflow (`.github/workflows/tag_release.yml`), which builds everything with the new version number and publishes it.
You can follow its progress on the [Actions](https://github.com/numbbo/coco-experiment/actions/workflows/tag_release.yml) tab on Github.


### Python

The workflow builds the Python source package and wheels for all platforms and Python versions in the matrix of `.github/workflows/test_python.yml` and tests the wheels.
It then publishes the source package and wheels to [PyPI](https://pypi.org/project/coco-experiment/) using [trusted publishing][trusted-publishing], so no PyPI token or manual upload is needed.
Publishing runs in the `release` environment; if that environment requires approval, approve the deployment on the workflow run page once the wheels have been built and tested.

Once the workflow has finished, check [PyPI](https://pypi.org/project/coco-experiment/), it should list the new release.

[trusted-publishing]: https://docs.pypi.org/trusted-publishers/


### Building additional wheels for an existing release

PyPI does not allow replacing files of a release, but it does allow adding new ones.
This lets you add wheels for new Python versions or platforms to a release after the fact.

1. If the new Python versions or platforms are not built for new releases yet, add them to the `python` (and, if needed, `buildplat`) matrix of the `python-wheels-build` job in `.github/workflows/test_python.yml`.
   Also add the new Python versions to the `python-wheels-test` and `python-wheels-example` jobs so the new wheels get tested.
   Merge the change into `main`.

1. Start the `tag-release` workflow manually for the release tag, either from the [Actions](https://github.com/numbbo/coco-experiment/actions/workflows/tag_release.yml) tab (*Run workflow*, branch `main`, tag e.g. `v2.8.2`) or with the [GitHub CLI](https://cli.github.com/):

   ```sh
   gh workflow run tag_release.yml -f tag=vX.Y.Z
   ```

   The workflow file is taken from `main`, the release itself is identified by the `tag` input.

1. If the `release` environment requires approval, approve the deployment once the wheels have been built and tested.

Instead of rebuilding the source package, the workflow downloads the sdist of that release from PyPI (checking its SHA256 digest) and builds all wheels in the matrix from it.
Before publishing, it drops every wheel PyPI already has for that version, so only the new wheels are uploaded to PyPI and attached to the existing GitHub release.
If nothing is new, nothing is uploaded.

Wheels count as the same if they share the Python and ABI tag and their platform tags only differ in the `manylinux`, `musllinux` or `macosx` version.
The release may have been built with an older `cibuildwheel` that used, e.g., `musllinux_1_1` instead of `musllinux_1_2`, and such a rebuilt wheel is not uploaded a second time.

### Other languages

The workflow attaches the source packages for C, Java, Matlab, Octave and Rust (`cocoex-<language>-X.Y.Z.zip`), together with the Python source package and wheels, to a [GitHub release](https://github.com/numbbo/coco-experiment/releases) for the tag.
