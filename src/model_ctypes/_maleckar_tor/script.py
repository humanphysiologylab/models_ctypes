import numpy as np
import pandas as pd
import ctypes
import matplotlib.pyplot as plt

import time
import os


dirname = '/home/andrey/WORK/Reentry/pypoptim/src/model_ctypes/_maleckar_tor/'
filename_so = os.path.join(dirname, 'maleckar.so')

filename_so_abs = os.path.abspath(filename_so)

model = ctypes.CDLL(filename_so_abs)

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


legend_constants = pd.read_csv(os.path.join(dirname, "legend_constants.csv"), index_col='name')['value']
legend_states = pd.read_csv(os.path.join(dirname, "legend_states.csv"), index_col='name')['value']

C = legend_constants.copy()
S = legend_states.copy()

n_beats = 10
t_sampling = 0.001
tol = 1e-6

chain_length = 50
v_threshold = 1e-1
t_safe = 1e-2

C['stim_period'] = 1
C['stim_duration'] = 0.001
C['stim_amplitude'] = -60

#C[7] = -80
stim_period = C['stim_period'] # C[5]
n_samples_per_stim = int(stim_period / t_sampling)

output = np.empty((n_samples_per_stim * n_beats + 1, len(S)))

status = model.run(S.values.copy(), C.values.copy(),
                   n_beats, t_sampling, tol, output)

print(status)

np.savetxt("./output.txt", output)
