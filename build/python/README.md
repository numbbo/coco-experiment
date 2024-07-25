<h1 align="center">
<img src="https://raw.githubusercontent.com/numbbo/coco-experiment/main/logo/coco-ex-300.webp" width="300">
</h1><br>

[![version](https://img.shields.io/pypi/v/coco-experiment?color=blue)](https://pypi.org/project/coco-experiment)
[![license](https://img.shields.io/pypi/l/coco-experiment)](https://pypi.org/project/coco-experiment)
[![build](https://img.shields.io/github/actions/workflow/status/numbbo/coco-experiment/build.yml?branch=main)][build]
[![download](https://img.shields.io/pypi/dm/coco-experiment.svg?label=PyPI%20downloads)][pypi]
[![paper](https://img.shields.io/badge/DOI-10.1080%2F10556788.2020.1808977-blue)][paper]

*coco-experiment* is a Python module implementing the experimental part of [COCO: A Platform for Comparing Continuous Optimizers in a Black-Box Setting][paper].

- **Website:** https://numbbo.github.io/coco/
- **Source code:** https://github.com/numbbo/coco-experiment
- **Bug reports:** https://github.com/numbbo/coco-experiment/issues

## Installation

The easiest way to install coco-experiment is using `pip` (preferably in a virtual environment):

```python
pip install -U coco-experiment
```

## Usage

Similar to scikit-learn, the name of the coco-experiment module is _not_ `coco-experiment` but rather `cocoex`.
To get started, take a look at the [simple example experiment][expy].

[build]: https://github.com/numbbo/coco-experiment/actions?query=branch%3Amain
[pypi]: https://pypi.org/project/coco-experiment/
[paper]: https://doi.org/10.1080/10556788.2020.1808977
[expy]: https://github.com/numbbo/coco-experiment/blob/main/build/python/example/example_experiment_complete.py
