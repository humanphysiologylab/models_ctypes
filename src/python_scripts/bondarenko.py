import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['figure.facecolor'] = 'white'
matplotlib.rcParams['savefig.dpi'] = 300
from matplotlib import cm

import seaborn as sns

from tqdm.auto import tqdm

import os
import ctypes
import gc

#  https://stackoverflow.com/a/37664693/13213091
def wrapped_ndptr(*args, **kwargs):
    base = np.ctypeslib.ndpointer(*args, **kwargs)
    def from_param(cls, obj):
        if obj is None:
            return obj
        return base.from_param(obj)
    return type(base.__name__, (base,), {'from_param': classmethod(from_param)})
DoubleArrayType_1D = wrapped_ndptr(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
DoubleArrayType_2D = wrapped_ndptr(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS')

#  https://stackoverflow.com/a/37664693/13213091
def wrapped_ndptr(*args, **kwargs):
    base = np.ctypeslib.ndpointer(*args, **kwargs)
    def from_param(cls, obj):
        if obj is None:
            return obj
        return base.from_param(obj)
    return type(base.__name__, (base,), {'from_param': classmethod(from_param)})
DoubleArrayType_1D = wrapped_ndptr(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
DoubleArrayType_2D = wrapped_ndptr(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS')

dirname = '../model_ctypes/_bondarenko'
filename_so = os.path.join(dirname, 'model.so')

filename_so_abs = os.path.abspath(filename_so)

model = ctypes.CDLL(filename_so_abs)

model.run.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), # double *S
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'), # double *C
    ctypes.c_int, # int n_beats
    ctypes.c_double, # double t_sampling
    ctypes.c_double, # double tol
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'), # double *output

    DoubleArrayType_2D, #np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'), # double *output_A
    DoubleArrayType_1D, #np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS') # double *output_t
    DoubleArrayType_1D, # double *stim_protocol
]

model.run.restype = ctypes.c_int


def run(S, C, n_beats, t_sampling, tol, output, *,
        output_A=None, t=None, stim_protocol=None):

    return model.run(S, C, n_beats, t_sampling, tol, output, output_A, t, stim_protocol)


legend_constants = pd.read_csv(os.path.join(dirname, "legend_constants.csv"), index_col='name')['value']
legend_states = pd.read_csv(os.path.join(dirname, "legend_states.csv"), index_col='name')['value']
# legend_algebraic = pd.read_csv(os.path.join(dirname, "legend_algebraic.csv"), index_col='name')
# legend_algebraic['value'] = 0.0
# legend_algebraic = legend_algebraic['value']

# model.initialize_states_default(legend_states.values, legend_constants.values)
# legend_states.to_csv(os.path.join(dirname, "legend_states.csv"))

S = legend_states.copy()
R = S * 0
C = legend_constants.copy()

t_sampling = 1

CL = 250
C['CL'] = CL
stim_period = C['CL']
n_samples_per_stim = int(stim_period / t_sampling)

n_beats = 1
tol = 1e-5

#%%time

output = np.zeros((n_samples_per_stim * n_beats + 1, len(S)))

Sc = S.values.copy()
Cc = C.values.copy()

#test_ctypes = ctypes.CDLL("./test.so")

status = run(Sc, Cc,
             n_beats, t_sampling, tol, output)

#print(status)

#output = pd.DataFrame(output, columns=legend_states.index)
