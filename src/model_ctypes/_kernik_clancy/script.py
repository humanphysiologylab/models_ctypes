import numpy as np
import pandas as pd
import ctypes
import matplotlib.pyplot as plt

import time
import os


filename_so = '/home/andrey/WORK/Reentry/pypoptim/src/model_ctypes/_koivumaki/koivumaki.so'
filename_so_abs = os.path.abspath(filename_so)

model = ctypes.CDLL(filename_so_abs)

model.initialize_states_default.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
]

model.initialize_states_default.restype = ctypes.c_void_p


model.initialize_constants_default.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
]

model.initialize_constants_default.restype = ctypes.c_void_p


model.run.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_int,
    ctypes.c_double,
    ctypes.c_double,
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS')
    #np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
    #np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
]

model.run.restype = ctypes.c_int


legend_constants = pd.read_csv("~/WORK/Reentry/pypoptim/src/model_ctypes/_koivumaki/legend_constants.csv")
legend_states = pd.read_csv("~/WORK/Reentry/pypoptim/src/model_ctypes/_koivumaki/legend_states.csv")

S = np.zeros(len(legend_states)) #  np.loadtxt("S.txt")
C = np.zeros(len(legend_constants)) #  np.loadtxt("C.txt")

model.initialize_states_default(S)
model.initialize_constants_default(C)
t_sampling = 0.001

C = pd.DataFrame(C, index=legend_constants['name'])
# C.loc['gCaL'] *= 5
C.loc['STIM_LEVEL'] = 1.0
C.loc['STIM_DURATION'] = 1e-3
C.loc['STIM_OFFSET'] = 0.0
C.loc['STIM_PERIOD'] = 1.

n_samples_per_stim = int(np.round(C.loc['STIM_PERIOD'] / t_sampling))

C = C.iloc[:, 0].values
#  print(S)
#  print(C)

n_beats = 2
tol = 1e-3

# chain_length = 50
# v_threshold = 1e-1
# t_safe = 1e-2

output = np.empty((n_samples_per_stim * n_beats + 1, len(S)))
# output_A = np.empty((n_samples_per_stim * n_beats + 1, 107))
# output_t = np.empty((n_samples_per_stim * n_beats + 1))

status = model.run(S.copy(), C, n_beats, t_sampling, tol, output)
print(status)

np.savetxt("output.txt", output)
