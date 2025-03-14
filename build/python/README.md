<h1 align="center">
    <table border="0">
  <td>
      <img src="https://raw.githubusercontent.com/numbbo/coco-experiment/main/logo/coco-ex-300.webp">          
  </td>
  <td>
      COmparing Continuous Optimisers COCO Experiment
</td>
</table>
</h1>

[![version](https://img.shields.io/pypi/v/coco-experiment?color=blue)](https://pypi.org/project/coco-experiment)
[![license](https://img.shields.io/pypi/l/coco-experiment)](https://pypi.org/project/coco-experiment)
[![build](https://img.shields.io/github/actions/workflow/status/numbbo/coco-experiment/build.yml?branch=main)][build]
[![download](https://img.shields.io/pypi/dm/coco-experiment.svg?label=PyPI%20downloads)][pypi]
[![paper](https://img.shields.io/badge/DOI-10.1080%2F10556788.2020.1808977-blue)][paper]

The `coco-experiment` Python module implements the experimental part of [COCO: A Platform for Comparing Continuous Optimizers in a Black-Box Setting][paper].
It provides several benchmark function suites and the data logging.

## Documentation

The main documentation pages for the `coco-experiment` package `cocoex` can be found at

- [getting-started](https://numbbo.it/getting-started#experiment)
- [API documentation](https://numbbo.github.io/coco-doc/apidocs/cocoex)
- [issue tracker and bug reports](https://github.com/numbbo/coco-experiment/issues)
- [source code page](https://github.com/numbbo/coco-experiment) on GitHub
- [COCO home page](https://coco-platform.org)

Below are installation instruction and some usage hints.

## Installation

The easiest way to install coco-experiment is using `pip` (preferably in a virtual environment):

```python
pip install -U coco-experiment
```

For the difficult way [look at this page][codepython] and start with `git clone https://github.com/numbbo/coco-experiment.git`.

## Usage

Similar to `scikit-learn`, _the module name_ (`coco-experiment`) **is not the import name**, hence usage starts with

```python
import cocoex
```

To get started, take a look at the [getting-started] page and at this [simple example experiment][expy]. To run an experiment in practice, take a look at the [complete example experiment][expy2].

[build]: https://github.com/numbbo/coco-experiment/actions?query=branch%3Amain
[pypi]: https://pypi.org/project/coco-experiment/
[paper]: https://doi.org/10.1080/10556788.2020.1808977
[getting-started]: https://coco-platform.org/getting-started/#experiment
[expy]: https://github.com/numbbo/coco-experiment/blob/main/build/python/example/example_experiment_simple.py
[expy2]: https://github.com/numbbo/coco-experiment/blob/main/build/python/example/example_experiment_complete.py
[codepython]: https://github.com/numbbo/coco-experiment/tree/main/build
