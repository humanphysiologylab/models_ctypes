import numpy as np
import pandas as pd
import itertools


def rastrigin(x, A=10):
    x = np.array(x)
    return sum(x**2 + A * (1 - np.cos(2 * np.pi * x)))


def uniform_vector(n=1):
    u = np.random.randn(n)
    return u / np.linalg.norm(u)


def cauchy_inverse_cdf(gamma):
    return gamma * np.tan(np.pi * (np.random.rand() - 0.5))


def calculate_RMSE(x, y):
    assert(len(x) == len(y))
    return np.sqrt(np.sum(((x - y)**2)) / len(x))


def update_array_from_kwargs(array, legend, **kwargs):
    """
    legend must contain `name` field with the same size as the array
    """

    series = pd.Series(data=array, copy=True,
                       index=legend['name'])

    for key in kwargs:
        if key in series:
            series[key] = kwargs[key]

    array_updated = series.values
    return array_updated


def get_value_by_key(array, legend, key):
    return dict(zip(legend['name'], array))[key]


def flatten_list(l):
    return list(itertools.chain(*l))


def batches_from_list(l, n_batches=1):
    return [l[i::n_batches] for i in range(n_batches)]


def value_from_bounds(bounds, log_scale=False):
    assert (len(bounds) == 2)
    assert (bounds[0] <= bounds[1])
    if log_scale:
        assert (bounds[0] > 0)
        bounds = np.log(bounds)
        value = np.exp(np.random.random() * np.ptp(bounds) + bounds[0])
    else:
        value = np.random.random() * np.ptp(bounds) + bounds[0]
    return value
