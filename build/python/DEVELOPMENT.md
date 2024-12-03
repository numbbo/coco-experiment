# Development Notes

Make sure that you
* install and activate the conda environment defined in `env.yaml` in the root directory of the repository.
  Alternatively, make sure you have all the require Python packages installed using some other method.
* run `scripts/fabricate` **before** you attempt to build or install the python package.

## Installing a development version

Just use `pip`.
Run
~~~sh
python -m pip install .
~~~
in the root directory of the cocoex Python package (`build/python` relativ to the repository root) to build a source distribution.

If you want to build with debug symbols, run
~~~sh
python -m pip install -e . \
    --no-build-isolation \
    -Csetup-args=-Dbuildtype=debug \
    -Cbuild-dir=build-dbg
~~~
instead.
This also enables `assert()` checking in the C code.

## Building a `sdist` archive

We use [meson-python](https://github.com/mesonbuild/meson-python) to build and install `cocoex`.
Unfortunately it does not support building source distributions if the package is not in the root of a git repository.
Therefore we cannot use it to build the sdist `tar.gz` archives.
Instead we use [hatchling](https://hatch.pypa.io/latest/) to build these.

Run
~~~sh
python -m hatchling build -t sdist -d dist/
~~~
in the root directory of the cocoex Python package (`build/python` relativ to the repository root) to build a source distribution.

## Building a wheel

You can use the [build](https://build.pypa.io/en/stable/)
Run
~~~sh
python -m build -w .
~~~
in the root directory of the cocoex Python package (`build/python` relativ to the repository root) to build a source distribution.

## Building binaries using meson

If you need to debug the meson build process, it might be useful to build the Python extensions directly.
Run
~~~sh
meson setup build
meson compile -C build
~~~
to setup and compile the extensions.

**Note**: It is not trivial to install this build using `meson install -C build` because it will try to install to the global Python `site-packages` directory!
