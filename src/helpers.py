import itertools
import ctypes
import os
import time

import numpy as np
import pandas as pd

from sklearn.metrics import mean_squared_error as MSE

def RMSE(x, y):
    return MSE(x, y, squared=False)


def strip_comments(code, comment_char='#'):
    lines = []
    for i, line in enumerate(code.split('\n')):
        items = line.split(comment_char)
        if len(items) >= 2:
            line = items[0]
        lines.append(line)
    code = '\n'.join(lines)
    return code


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
    return np.sqrt(np.mean(((x - y)**2)))


def calculate_RMSE_balanced(x, y):
    assert(len(x) == len(y))
    x = (x - y.min(axis=0)) / y.ptp(axis=0)
    y = (y - y.min(axis=0)) / y.ptp(axis=0)
    return np.sqrt(np.sum(((x - y)**2)) / len(x))


def autoscaling(signal_to_scale, signal_reference):
    def scalar_multiplications(a, b):
        assert (len(a) == len(b))
        assert (a.size == b.size)
        coefficients = np.array([np.dot(a, b), np.sum(a), np.sum(b), np.sum(a ** 2), len(a)])
        return coefficients

    c = scalar_multiplications(signal_to_scale, signal_reference)

    if c[1] == 0 or c[1] * c[1] - c[4] * c[3] == 0:
        alpha = 0
        beta = 0
    else:
        beta = (c[0] * c[1] - c[2] * c[3]) / (c[1] * c[1] - c[4] * c[3])
        alpha = (c[2] - beta * c[4]) / c[1]

    signal_scaled = signal_to_scale * alpha + beta
    rmse = calculate_RMSE(signal_scaled, signal_reference)

    return signal_scaled, rmse, (alpha, beta)


def calculate_composite_RMSE_V_CaT(x, y):
    # x -- model, y -- experiment
    assert (len(x) == len(y))

    v_model, cat_model = x.T
    v_exp, cat_exp = y.T

    cat_model = (cat_model - cat_model.min(axis=0)) / cat_model.ptp(axis=0)  # to balance V and CaT

    rmse_v = calculate_RMSE_balanced(v_model, v_exp)    # v_exp --> [0, 1]
    cat_exp_scaled, rmse_cat, coeffs = autoscaling(signal_to_scale=cat_exp,
                                                   signal_reference=cat_model)

    rmse_total = rmse_v + rmse_cat

    return rmse_total


# # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # #
def ma(x, N):
    return np.convolve(x, np.ones(N)/N, mode='valid')

def calculate_mean_abs_noise(array, N=31):
    array_ma = np.apply_along_axis(func1d=ma, axis=0, arr=array, N=N)
    array_valid = array[N//2: N//2 + len(array_ma)]
    noise = np.mean(np.abs(array_valid - array_ma), axis=0)
    return noise

def calculate_RMSE_weightened(x, y, weights):
    return np.sum(MSE(x, y, squared=False, multioutput='raw_values') * weights)
# # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # #


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


def flatten_iterable(x):
    return list(itertools.chain(*x))


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


def argmax_list_of_dicts(l, key):
    return max(enumerate(l), key=lambda x: x[1][key])[0]


def find_index_first(seq, condition):
    #  https://stackoverflow.com/a/8534381/13213091
    return next((i for i, x in enumerate(seq) if condition(x)), None)


def create_model_from_so(filename_so):

    #filename_so = '_maleckar/maleckar.so'
    #dirname_current_abs = os.path.abspath(os.path.dirname(__file__))
    #filename_so_abs = os.path.join(dirname_current_abs, filename_so)
    #model = ctypes.CDLL(filename_so_abs)
    model = ctypes.CDLL(filename_so)

    model.initConsts.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
    ]

    model.initConsts.restype = ctypes.c_void_p


    model.run.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    ]

    model.run.restype = ctypes.c_int


    model.run_chain.argtypes = [
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_int,
        ctypes.c_double,
        ctypes.c_double,
        np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    ]

    model.run_chain.restype = ctypes.c_int

    return model


class Timer:

    def __init__(self):
        self.times = dict()
        self.total = None

    def start(self, name):
        self.times[name] = -time.time()

    def end(self, name):
        self.times[name] += time.time()

    def report(self, sort=False):
        total = sum(self.times.values())
        self.times['total'] = total
        times = dict(sorted(self.times.items(), key=lambda x: x[1], reverse=True) if sort else self.times)
        s = "\n".join(f"{name}: {times[name]:.6f} {100 * times[name] / total:.2f}%" for name in times)
        del self.times['total']
        return s

    def clear(self):
        del self.times
        self.times = dict()
