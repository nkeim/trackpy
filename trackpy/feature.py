from __future__ import division
import numpy as np
import math

import numba
from numba import double, int_, bool_
# Both the explicit signature below, and the delayed compile,
# produce the slow import.
#@numba.jit(double[:,:](int_[:,:], int_[:,:], int_, double[:,:],
#                int_, bool_, int_[:], int_[:,:], int_[:,:],
#                double[:,:], double[:,:], int_),)
@numba.jit
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

            # If we're off by more than half a pixel in any direction, move.
            do_move = False
            if allow_moves:
                for dim in range(2):
                    if off_center[dim] > SHIFT_THRESH:
                        do_move = True

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


