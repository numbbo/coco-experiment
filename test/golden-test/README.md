# Golden master tests for the core COCO framework

These tests record the observable behaviour of COCO and check that it does not change.
They are meant as a safety net for refactorings that should not change any results, such as the port to C++.

The driver `golden.c` only uses the public C interface from `coco.h`.

- **values** scenarios evaluate every problem of a suite (default instances) at its initial solution, three random points in the region of interest, one point that is likely outside of it, the initial solution with one coordinate set to NaN, +Inf and -Inf, and the optimal solution.
  They record the problem metadata, objective and constraint values and the evaluation counters.
  The output has one line per (function, dimension) pair with a hash over all values of its instances.
- **observe** scenarios run a deterministic random search (followed by a NaN and an Inf solution) with an observer attached and record a hash of every file written by the logger.
  Only `coco_version` is masked because it changes with every commit.

All comparisons are exact.

## Running the tests

1. Fabricate the sources: `python ../../scripts/fabricate`
1. Configure the build: `cmake -B build`
1. Build: `cmake --build build`
1. Run the tests: `ctest --test-dir build`

## Platforms

Results depend on the platform's math library, so references are kept per platform in `reference/<system>-<processor>`.
If there are no references for your platform, the tests are disabled and CMake prints a warning.
To create or update the references for your platform, run

~~~
cmake --build build --target golden-update
~~~

and then configure again so that the tests are registered.
Review the changes to `reference/` before you commit them; any change means that COCO's results changed.

The library is always compiled with `-O2` because the optimizer replaces some `libm` calls (e.g. `pow(10, x)` by `exp10(x)`), which can change results in the last bit.

## Finding out what changed

When a test fails, it prints the lines that differ.
To see which value changed, dump all values with both the old and the new version and compare the dumps:

~~~
./build/golden values bbob --dump > new.txt
diff old.txt new.txt
~~~
