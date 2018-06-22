from __future__ import print_function
import os
import inspect

# (copied from the BEAST config.py file)


# Set to use some C code instead of pure python to speed up the computations.
# If False, only numpy and python code are used.
#__WITH_C_LIBS__ = True
__WITH_C_LIBS__ = False

#numexpr -- optimized multi-threaded numpy evaluations
__USE_NUMEXPR__ = True

# Default number of threads to use when parallel computing
#__NTHREADS__ = 6
__NTHREADS__ = 1

#directories
__ROOT__ = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])
