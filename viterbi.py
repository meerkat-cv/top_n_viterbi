import ctypes
import numpy as np

lib = ctypes.CDLL('./build/libviterbi.so')

print_numpy_float32_array = lib.print_numpy_float32_array
print_numpy_float32_array.argtypes = [
    np.ctypeslib.ndpointer(ctypes.c_float, flags='C_CONTIGUOUS'),
    ctypes.c_int,
]

# ocr_viterbi_topk = lib.ocr_viterbi_topk


test_array = np.asarray([1.2, 4.5, 2.4, -89.4], dtype=np.float32)
print_numpy_float32_array(test_array, test_array.shape[0])
print('new array', test_array)
