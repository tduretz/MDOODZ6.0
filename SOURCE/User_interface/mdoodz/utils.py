#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 31 10:52:49 2019

@author: abauville
"""

use_numba = True
try:
    from numba import jit #, prange
    use_numba = True
except ImportError:
    print("Warning: the package numba was not found. just-in-time compilation is off")
    use_numba = False
    
def maybe_numba(use_numba):
    return jit(nopython=True) if use_numba else lambda x:x



