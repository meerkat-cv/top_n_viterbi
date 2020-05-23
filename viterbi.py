import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer

lib = ctypes.CDLL('./build/libviterbi.so')

ocr_viterbi_topk = lib.ocr_viterbi_topk
ocr_viterbi_topk.argtypes = [
    ndpointer(ctypes.c_float, flags='C_CONTIGUOUS'),
    ndpointer(ctypes.c_float, flags='C_CONTIGUOUS'),
    ndpointer(ctypes.c_float, flags='C_CONTIGUOUS'),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_int,
    ndpointer(ctypes.c_int, flags='C_CONTIGUOUS'),
]


pi = np.asarray([0.5,0.5], dtype=np.float32)
a = np.asarray([[0.5, 0.5], [0.5, 0.5]], dtype=np.float32)
b = np.asarray([[0.2, 0.4, 0.8], [0.8, 0.6, 0.2]], dtype=np.float32)
top_n = 4

paths = np.zeros((top_n, b.shape[1]), dtype=np.int32)
ocr_viterbi_topk(pi, a, b, b.shape[0], b.shape[1], top_n, paths)

print('paths', paths)

lib.free_variables()
