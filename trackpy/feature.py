# Copyright 2012 Daniel B. Allan
# dallan@pha.jhu.edu, daniel.b.allan@gmail.com
# http://pha.jhu.edu/~dallan
# http://www.danallan.com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses>.


from __future__ import division
#import warnings
import numpy as np
#import pandas as pd
#from scipy import ndimage
#from scipy import stats
#from pandas import DataFrame, Series
#import matplotlib.pyplot as plt  # for walkthrough

def _get_numba_refine_locals():
    """Establish types of local variables in _numba_refine(), in a way that's safe if there's no numba."""
    try:
        from numba import double, int_, bool_
    except ImportError:
        return {}
    else:
        return dict(SHIFT_THRESH=double, GOOD_ENOUGH_THRESH=double,
                    square0=int_, square1=int_, square_size=int_, Rg_=double,
                    ecc_=double, mass_=double, ecc1=double, ecc2=double,
                    eccsq=double,
                    signal_=double, allow_moves=bool_, do_move=bool_,
                    center_px=int_, raw_px=int_, px=int_,
                    oc=double, result=double[:,:],
                    final_coords=double[:,:], mass=double[:],
                    Rg=double[:], ecc=double[:], signal=double[:],
                    coord=double[:], cm_n=double[:],
                    cm_i=double[:], off_center=double[:],
                    new_coord=int_[:], )

import numba
import math
from numba import double, int_, bool_
numba.config.DEBUG = 1
numba.config.OPT = 0
#@try_numba_autojit(locals=_get_numba_refine_locals())
#@try_numba_autojit
@numba.jit(double[:,:](int_[:,:], int_[:,:], int_, double[:,:],
                int_, bool_, int_[:], int_[:,:], int_[:,:],
                double[:,:], double[:,:], int_),
           locals=_get_numba_refine_locals())
def _numba_refine(raw_image, image, radius, coords, max_iterations,
                  characterize, shape, mask, r2_mask, cmask, smask, N):
    SHIFT_THRESH = 0.6
    GOOD_ENOUGH_THRESH = 0.01

    square_size = 2*radius + 1

    # Declare arrays that we will fill iteratively through loop.
    #N = coords.shape[0]
    final_coords = np.zeros((N, 2), dtype=np.double)
    mass = np.zeros(N, dtype=np.double)
    Rg = np.zeros(N, dtype=np.double)
    ecc = np.zeros(N, dtype=np.double)
    signal = np.zeros(N, dtype=np.double)
    coord = np.zeros(2, dtype=np.double)

    # Buffer arrays
    cm_n = np.zeros(2, dtype=np.double)
    cm_i = np.zeros(2, dtype=np.double)
    off_center = np.zeros(2, dtype=np.double)
    new_coord = np.zeros(2, dtype=np.int64)

    # "Declare" registers
    square0 = 0
    square1 = 0
    Rg_ = 0.
    ecc_ = 0.
    mass_ = 0.
    ecc1 = 0.
    ecc2 = 0.
    signal_ = 0.
    allow_moves = True
    do_move = False
    center_px = image[0,0]
    raw_px = raw_image[0,0]
    px = image[0,0]
    oc = 0.

    for feat in range(N):
        # Define the circular neighborhood of (x, y).
        for dim in range(2):
            coord[dim] = coords[feat, dim]
            cm_n[dim] = 0.
        square0 = int(coord[0] - radius)
        square1 = int(coord[1] - radius)
        mass_ = 0.0
        for i in range(square_size):
            for j in range(square_size):
                if mask[i, j] != 0:
                    px = image[square0 + i, square1 + j]
                    cm_n[0] += px*i
                    cm_n[1] += px*j
                    mass_ += px

        for dim in range(2):
            cm_n[dim] /= mass_
            cm_i[dim] = cm_n[dim] - radius + coord[dim]
        allow_moves = True
        for iteration in range(max_iterations):
            for dim in range(2):
                off_center[dim] = cm_n[dim] - radius
            for dim in range(2):
                if off_center[dim] > GOOD_ENOUGH_THRESH:
                    break  # Proceed through iteration. # FIXME
                break  # Stop iterations. # FIXME

            # If we're off by more than half a pixel in any direction, move.
            do_move = False
            if allow_moves:
                for dim in range(2):
                    if off_center[dim] > SHIFT_THRESH:
                        do_move = True
            #do_move = True # FIXME: Why do we always need to move?

            if do_move:
                # In here, coord is an integer.
                for dim in range(2):
                    new_coord[dim] = int(round(coord[dim]))
                    oc = off_center[dim]
                    if oc > SHIFT_THRESH:
                        new_coord[dim] += 1
                    elif oc < - SHIFT_THRESH:
                        new_coord[dim] += -1
                    # Don't move outside the image!
                    if new_coord[dim] < radius:
                        new_coord[dim] = radius
                    upper_bound = shape[dim] - radius - 1
                    if new_coord[dim] > upper_bound:
                        new_coord[dim] = upper_bound
                # Update slice to shifted position.
                square0 = new_coord[0] - radius
                square1 = new_coord[1] - radius
                for dim in range(2):
                     cm_n[dim] = 0.

            # If we're off by less than half a pixel, interpolate.
            else:
                break
                # TODO Implement this for numba.
                # Remember to zero cm_n somewhere in here.
                # Here, coord is a float. We are off the grid.
                # neighborhood = ndimage.shift(neighborhood, -off_center,
                #                              order=2, mode='constant', cval=0)
                # new_coord = np.float_(coord) + off_center
                # Disallow any whole-pixels moves on future iterations.
                # allow_moves = False

            # cm_n was re-zeroed above in an unrelated loop
            mass_ = 0.
            for i in range(square_size):
                for j in range(square_size):
                    if mask[i, j] != 0:
                        px = image[square0 + i, square1 + j]
                        cm_n[0] += px*i
                        cm_n[1] += px*j
                        mass_ += px

            for dim in range(2):
                cm_n[dim] /= mass_
                cm_i[dim] = cm_n[dim] - radius + coord[dim]
                coord[dim] = new_coord[dim]
        # matplotlib and ndimage have opposite conventions for xy <-> yx.
        final_coords[feat, 0] = cm_i[1]
        final_coords[feat, 1] = cm_i[0]

        # Characterize the neighborhood of our final centroid.
        mass_ = 0.
        Rg_ = 0.
        ecc1 = 0.
        ecc2 = 0.
        signal_ = 0.
        for i in range(square_size):
            for j in range(square_size):
                if mask[i, j] != 0:
                    px = image[square0 + i, square1 + j]
                    mass_ += px
                    # Will short-circuiting if characterize=False slow it down?
                    if not characterize:
                        continue
                    else:
                        Rg_ += r2_mask[i, j]*px
                        ecc1 += cmask[i, j]*px
                        ecc2 += smask[i, j]*px
                        raw_px = raw_image[square0 + i, square1 + j]
                        if raw_px > signal_:
                            signal_ = px
        Rg_ = math.sqrt(Rg_/mass_)
        mass[feat] = mass_
        if characterize:
            Rg[feat] = Rg_
            center_px = image[square0 + radius, square1 + radius]
            eccsq = ecc1**2 + ecc2**2
            ecc_ = math.sqrt(eccsq)/(mass_ - center_px + 1.0e-6)
        ecc[feat] = ecc_
        signal[feat] = signal_  # black_level subtracted later

    if not characterize:
        result = np.column_stack([final_coords, mass])
    else:
        result = np.column_stack([final_coords, mass, Rg, ecc, signal])
    return result


