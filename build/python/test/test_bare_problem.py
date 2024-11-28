import pytest
import numpy as np

from cocoex import BareProblem
from cocoex.exceptions import NoSuchSuiteException


def test_bad_suite():
    with pytest.raises(NoSuchSuiteException):
        BareProblem("bad_suite", 1, 1, 1)


def test_bbob():
    for fid in range(1, 25):
        for dimension in [2, 3, 4, 5, 7, 11, 13, 17, 18, 19, 20, 31, 40]:
            x0 = np.zeros(dimension)
            for instance in [0, 1, 10, 100, 1000]:
                fn = BareProblem("bbob", fid, dimension, instance)
                assert str(fn) == f"bbob_f{fid:03d}_i{instance:02d}_d{dimension:02d}"
                assert fn(x0) >= fn.best_value()
                assert not fn.is_tainted
                assert np.isclose(fn.best_value(), fn(fn.best_parameter()))
                assert fn.is_tainted


def test_list():
    fn = BareProblem("bbob", 1, 4, 1)
    assert fn([1, 2, 3, 4]) >= fn([1.0, 2.0, 3.0, 4.0])


def test_multiple_parameters():
    for n in [1, 2, 10, 100, 200]:
        fn = BareProblem("bbob", 1, 4, 1)
        X = np.random.uniform(-5, 5, size=(n, 4))
        y = fn(X)
        assert len(y) == n
        assert np.all(y >= fn.best_value())
