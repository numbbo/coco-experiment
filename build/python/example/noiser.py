import collections
import warnings
import math
import json
import numpy as np  # for checking finite values

_frozen_noise = True
math_floor = int
parameters_file = 'noiser_parameters.json'

def rand(x, i=0):
    """return a uniform random value seeded with ``x[:2]`` and `i`"""
    return _unif(1, _seed(x, i))[0]
def randn(x, i=0):
    """return a normal random value seeded with ``x[:2]`` and `i`"""
    return _gauss(1, _seed(x, i))[0]
def randc(x, i=0):
    """return a Cauchy random value seeded with ``x[:2]`` and `i`.

    We observe that median |Z| = 1 and p(pi x Z / 2 > a) = 1/a.
    """
    return _cauchy(1, _seed(x, i))[0]

class NoisifyProblem:
    """noise wrapper for an observed problem.

    The noise is applied with different probabilities for

    - subtracting a value (leading to good outliers)
    - adding a value (leading to bad outliers)
    - adding a small `epsilon` normally distributed noise

    The noise is frozen in x and y: whether it's applied and its
    sampled value depend deterministically on the input x and i.

    The outlier noise distribution is a Cauchy distribution multiplied by
    pi/2, which generates values larger than k with a probability awfully
    close to k**-1 when k >= 10.

    Example::

        import cocoex  # experimentation module
        import noiser  # noise.py needs to be in the Python path
        import scipy  # to define the solver to be benchmarked

        ### input
        suite_name = "bbob"
        fmin = scipy.optimize.fmin  # optimizer to be benchmarked

        ### prepare
        suite = cocoex.Suite(suite_name, "", "")
        observer = cocoex.Observer(suite_name, "")

        ### go
        for problem in suite:  # this loop may take several minutes or more
            problem.observe_with(observer)  # generates the data for cocopp
            problem = noiser.NoisifyProblem(problem)
            fmin(problem, problem.initial_solution, disp=False)

    Details: the random number generators can be passed as argument and
    need to obey the interface of `noiser.rand`. ``rands[0]`` is assumed to
    be uniform in [0,1], ``rands[1]`` is used to sample the Gaussian noise,
    ``abs(rands[2])`` is used to sample added noise and ``abs(rands[-1])``
    is used to sample subtracted noise.

    See https://github.com/numbbo/coco-experiment/blob/main/build/python/example/example_experiment_complete.py

    """
    def __init__(self, problem, p_add=0.2, p_subtract=0.0, p_epsilon=0, epsilon=1e-4,
                 rands=(rand, randn, randc)):
        self._params = {k: v for k, v in locals().items() if k != 'self'}
        self._problem = problem
        if sum(self._params[p] for p in ('p_subtract', 'p_add')) > 1:
            warnings.warn("p_subtract={0} + p_add={1} > 1, hence p_add is interpreted"
                          " as ``1-p_subtract``"
                          .format(self._params['p_subtract'], self._params['p_add']))
    @property
    def parameters(self):
        return {k: v for (k, v) in self._params.items() if k.startswith(('p_', 'eps'))}

    def _dump_params(self):
        with open(parameters_file, 'wt') as fp:
            json.dump(self.parameters, fp)

    def _load_params(self, filename=parameters_file, warn_errors=(FileNotFoundError,), ignore_errors=()):
        try:
            with open(parameters_file, 'rt') as fp:
                return json.load(fp)
        except warn_errors as e:
            warnings.warn('NoisifyProblem._load_params raised the exception "{0}"'.format(e))
        except ignore_errors:
            pass
        return {}

    def _set_params(self, filename=parameters_file, **kwargs):
        """set parameters from file `filename`"""
        self._params.update(self._load_params(filename), **kwargs)

    def __getattr__(self, name):
        """pretend to be a COCO `problem` as passed in `__init__`"""
        return getattr(self._problem, name)

    def __call__(self, x):
        return self._fnoise(x) + self._problem(x)

    def constraint(self, x):
        constraint = self._problem.constraint(x)
        for i, c in enumerate(constraint):
            constraint[i] += self._cnoise([x[i % 2], c])
        return constraint

    def rand1(self, x, fac=1):
        return self._params['rands'][0](x, fac)
    def rand2(self, x, fac=1):
        return self._params['rands'][1](x, fac)
    def rand3(self, x, fac=1):
        return self._params['rands'][2](x, fac)

    def _cnoise(self, x):
        return self._fnoise(x)  # for the time being

    def _fnoise(self, x):
        """additive x-dependent noise"""
        n = 0  # no noise by default
        if self._params['p_epsilon'] > 0 and self.rand1(x, 3) < self._params['p_epsilon']:
            n = self._params['epsilon'] * self.rand2(x)
        assert np.isfinite(n)
        if all(self._params[p] <= 0 for p in ('p_add', 'p_subtract')):
            return n
        r = self.rand1(x, 2)  # a uniform random number
        assert 0 <= r <= 1, r
        i = 4
        for s, p, rand in [(1, 'p_add', 2),
                           (-1, 'p_subtract', -1)]:
            p = self._params[p]
            if r < p:
                assert p > 0, (r, p)
                nn = self._params['rands'][rand](x, i)
                if nn < 0:
                    nn *= -1
                return n + s * 3.141592653589793 / 2 * nn
                        # self.randn(x, 2) / self.randn(x, 12)
                # was: return n + s * 2 * r / self._params[p] / self.rand1(x)
            r = 1 - r  # prevent sampling another r
            i += 1
        assert np.isfinite(n)
        return n

_seed_weights = [1.23468, 2.34579]  # otherwise, coordinate search would be noisefree
_seeds = collections.deque(maxlen=99)  # for the record only
def _seed(x, i=0):
    """return a scalar seed using the first two entries of vector `x`"""
    if not _frozen_noise:
        import random
        return 1e4 / (random.random() + 1e-22)
    freezer = (i + 1) * (1 + _seed_weights[0] * float(x[0]) + _seed_weights[1] * float(x[1]))
    if freezer < 0:
        freezer = -freezer
    # s = 1 + int(repr(freezer).replace('.', '').split('e')[0][-1:-15:-1])
    s = freezer + 1 / (freezer + 136247e-22)
    while s > 1e21:
        s /= 9
    _seeds.appendleft(s)
    if not np.isfinite(s):
        warnings.warn("seed({0}) = {1} is not finite. Using 1 instead"
                      .format(x, s))
        return 1
    return s
def _unif(N, inseed):
    """Generates N uniform numbers with starting seed inseed.

    This is copy-paste from the bbobbenchmarks Python implementation with
    minor adaptation.
    """

    # initialization
    if inseed < 0:  # np.abs took half the overall time without showing up in the profiler
        inseed = -inseed
    if inseed < 1.:
        inseed += 1.

    rgrand = 32 * [0.]
    aktseed = inseed
    for i in range(39, -1, -1):
        tmp = math_floor(aktseed / 127773.)
        aktseed = 16807. * (aktseed - tmp * 127773.) - 2836. * tmp
        if aktseed < 0:
            aktseed = aktseed + 2147483647.
        if i < 32:
            rgrand[i] = aktseed
    aktrand = rgrand[0]

    # sample numbers
    r = N * [0.]
    N0 = 2  # the first is visibly biased to be < 0.05
    for i in range(-N0, N):
        tmp = math_floor(aktseed / 127773.)
        aktseed = 16807. * (aktseed - tmp * 127773.) - 2836. * tmp
        if aktseed < 0:
            aktseed = aktseed + 2147483647.
        tmp = math_floor(aktrand / 67108865.)
        aktrand = rgrand[tmp]
        rgrand[tmp] = aktseed
        if i >= 0:
            r[i] = aktrand / 2.147483647e9
    # r = np.asarray(r)
    if not all(r):
        warnings.warn('zero sampled {0} times, set to 1e-99'.format(sum(ri == 0 for ri in r)))
        for i in range(len(r)):
            if r[i] == 0:
                r[i] = 1e-99
    return r
def _gauss(N, seed):
    """Samples N standard normally distributed numbers
    being the same for a given seed

    """
    r = _unif(2 * N, seed)
    g = [math.sqrt(-2 * math.log(r[i])) * math.cos(2 * math.pi * r[i+N]) for i in range(N)]
    if not all(g):
        for i in range(len(g)):
            if g[i] == 0:
                g[i] = 1e-99
    return g
def _cauchy(N, seed):
    """Samples N standard normally distributed numbers
    being the same for a given seed

    """
    r = _unif(4 * N, seed)
    g = [math.sqrt(-2 * math.log(r[i])) * math.cos(2 * math.pi * r[i+N]) for i in range(N)]
    return [g[i] / max((1e-21,
            math.fabs(math.sqrt(-2 * math.log(r[i+2*N])) * math.cos(2 * math.pi * r[i+3*N]))))
            for i in range(N)]
# def _gauss(N, seed):  # with numpy
#     """Samples N standard normally distributed numbers
#     being the same for a given seed

#     """
#     r = np.asarray(_unif(2 * N, seed))
#     g = np.sqrt(-2 * np.log(r[:N])) * np.cos(2 * np.pi * r[N:2*N])
#     if not np.all(g):
#         g[g == 0] = 1e-99
#     return g
