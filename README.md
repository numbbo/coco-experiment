<h1 align="center">
    <table border="0">
  <td>
    <img src="https://raw.githubusercontent.com/numbbo/coco-experiment/main/logo/coco-ex-300.webp" width="300">
  </td>
  <td>
      <A href="https://coco-platform.org">COmparing Continuous Optimisers</A> (COCO) <A href="https://numbbo.github.io/coco-doc/apidocs/cocoex">experiment</A>
  </td>
</table>
</h1>

[![build](https://img.shields.io/github/actions/workflow/status/numbbo/coco-experiment/build.yml?branch=main)][buildaction]
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2594848.svg)][paper]
[![PyPI Downloads](https://static.pepy.tech/badge/coco-experiment/month)](https://pepy.tech/projects/coco-experiment)

> [!IMPORTANT]
>
> This repository contains the source files for the experiment part of the `coco` framework.
> If you don't want to *extend* the framework, you probably don't need this!
> See instead the new [COCO documentation](https://numbbo.it/getting-started) and use the language bindings of your choice from the package repository for your language (e.g. [PyPI](https://pypi.org/project/coco-experiment) for Python, [crates.io](https://crates.io) for Rust, ..., or otherwise see the [build folder](https://github.com/numbbo/coco-experiment/tree/main/build)).

This repository implements the experiment part of [COCO: A Platform for Comparing Continuous Optimizers in a Black-Box Setting][paper], comparing not only continuous optimisers.
This code generates the data used within the [COCO framework](https://coco-platform.org).

The documentation for COCO experiment can be found at

- [getting-started](https://numbbo.it/getting-started#experiment)
- [C API](https://numbbo.github.io/coco-doc/C) documentation
- [Python API](https://numbbo.github.io/coco-doc/apidocs/cocoex) documentation
- [Issue tracker and bug reports](https://github.com/numbbo/coco-experiment/issues)
- [source code page](https://github.com/numbbo/coco-experiment) on GitHub
- [COCO home page](https://coco-platform.org)


[buildaction]: https://github.com/numbbo/coco-experiment/actions?query=branch%3Amain
[paper]: https://doi.org/10.5281/zenodo.2594848
