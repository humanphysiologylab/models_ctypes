import numpy as np
import ctypes
import matplotlib.pyplot as plt

import time
import os


filename_so = '/home/andrey/WORK/Reentry/pypoptim/src/model_ctypes/_maleckar/maleckar.so'
filename_so_abs = os.path.abspath(filename_so)

maleckar = ctypes.CDLL(filename_so_abs)

maleckar.initConsts.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS')
]

maleckar.initConsts.restype = ctypes.c_void_p


maleckar.run.argtypes = [
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.c_int,
    ctypes.c_double,
    ctypes.c_double,
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='C_CONTIGUOUS'),
]

maleckar.run.restype = ctypes.c_int


maleckar.run_chain.argtypes = [
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

maleckar.run_chain.restype = ctypes.c_int


S = np.zeros(30) #  np.loadtxt("S.txt")
C = np.zeros(51) #  np.loadtxt("C.txt")

maleckar.initConsts(C, S)

n_beats = 10
t_sampling = 0.001
tol = 1e-3

chain_length = 50
v_threshold = 1e-1
t_safe = 1e-2

C[7] = -80
stim_period = C[5]
n_samples_per_stim = int(stim_period / t_sampling)

output = np.empty((n_samples_per_stim * n_beats + 1, len(S)))

status = maleckar.run_chain(S.copy(), C, chain_length, v_threshold, t_safe, n_beats, t_sampling, tol, output)

np.savetxt("./output.txt", output)
