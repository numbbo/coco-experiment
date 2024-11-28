# -*- mode: cython -*-
#cython: language_level=3, boundscheck=False, c_string_type=str, c_string_encoding=ascii

import numpy as np
cimport numpy as np

from .exceptions import NoSuchProblemException, NoSuchSuiteException

np.import_array()


cdef extern from 'coco.h':
    ctypedef struct coco_problem_t:
        pass

    char *coco_problem_get_id(coco_problem_t *p)
    void coco_evaluate_function(coco_problem_t *p, double *x, double *y)
    void coco_problem_free(coco_problem_t *p)
    void bbob_problem_best_parameter_print(coco_problem_t *p)
    
    int coco_problem_get_best_parameter(const coco_problem_t *problem, double *best_parameter)
    int coco_problem_is_tainted(const coco_problem_t *problem)

    int coco_problem_get_best_parameter(const coco_problem_t *problem, double *best_parameter)

# IMPORTANT: These functions are *not* declared public in coco.h so we have to
# explicitly declare it as an external function. Otherwise Cython will *not*
# add a declaration to the generated source files.
cdef extern coco_problem_t *coco_get_bbob_problem(size_t function, size_t
                                                  dimension, size_t instance)
cdef extern double coco_problem_get_best_value(coco_problem_t *p)


cdef class BareProblem:
    """A bare benchmark function from one of the available suites.

    Examples
    --------

    >>> import numpy as np

    Create a 13 dimensional sphere function

    >>> fn = BareProblem("bbob", 1, 13, 1)
    >>> fn
    BareProblem('bbob', 1, 13, 1)

    We can also get a short mnemonic name for the function

    >>> str(fn)
    'bbob_f001_i01_d13'

    And of course evaluate it

    >>> x = np.zeros(fn.dimension)
    >>> fn(x)
    124.61122368000001
    """

    cdef coco_problem_t *_problem
    cdef readonly char* suite_name
    cdef readonly char* id
    cdef readonly int function
    cdef readonly int dimension
    cdef readonly int instance

    def __init__(self, suite_name: str, function: int, dimension: int, instance: int):
        """
        Create a bare benchmark problem from one of the COCO suites.

        Parameters
        ----------
        suite_name
            Name of benchmark suite ("bbob" only currently)

        function
            ID of function from suite

        dimension
            Dimension of instantiated function

        instance
            Instance ID of instantiated function


        Raises
        ------
        NoSuchSuiteException
           If the `suite_name` is not known or not yet supported

        NoSuchProblemException
          If no problem with the given `function`, `dimension` and `instance` exists in
          the given `suite_name`.
        """
        self.suite_name = suite_name
        self.function = function
        self.dimension = dimension
        self.instance = instance
        self._problem = NULL
        if suite_name == "bbob":
            self._problem = coco_get_bbob_problem(function, dimension, instance)
            if self._problem == NULL:
                # FIXME: Possibly extend Exception to include dimension and instance?
                raise NoSuchProblemException(suite_name, function)
        else:
            raise NoSuchSuiteException(suite_name)

        self.id = coco_problem_get_id(self._problem)

    def __del__(self):
        if self._problem != NULL:
            coco_problem_free(self._problem)

    def best_value(self):
        """Return the best (lowest) possible function value"""
        return coco_problem_get_best_value(self._problem)

    def best_parameter(self):            
        cdef double[::1] x
        x = np.zeros(self.dimension, dtype=np.float64)

        res = coco_problem_get_best_parameter(self._problem, &x[0])
        if res == 0: # Failure
            return None # FIXME: Raise an exception?
        else:
            return np.asarray(x)

    @property
    def is_tainted(self):
        is_tainted = coco_problem_is_tainted(self._problem)
        return is_tainted != 0
    
    def __str__(self):
        return self.id

    def __repr__(self):
        return f"BareProblem('{self.suite_name}', {self.function}, {self.dimension}, {self.instance})"

    def __call__(self, x):
        cdef double[::1] xi
        cdef double[:, ::1] X
        cdef Py_ssize_t N, D, i
        cdef double[::1] y_view
        cdef double y
        if isinstance(x, list):
            x = np.array(x, dtype=float)
        if x.ndim == 1:
            # Evaluate a single parameter
            xi = np.array(x, dtype=float)
            coco_evaluate_function(self._problem, &xi[0], &y)
            return y
        elif x.ndim == 2:
            # Evaluate several parameters at once
            X = x
            N = X.shape[0]
            D = X.shape[1]
            Y = np.zeros(N, dtype=np.float64)
            y_view = Y
            for i in range(N):
                xi = X[i, :]
                coco_evaluate_function(self._problem, &xi[0], &y)
                y_view[i] = y
            return Y
        else:
          return None

