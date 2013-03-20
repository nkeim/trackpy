#Copyright 2013 Thomas A Caswell
#tcaswell@uchicago.edu
#http://jfi.uchicago.edu/~tcaswell
#
#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or (at
#your option) any later version.
#
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses>.
from __future__ import division

import numpy as np
import numpy.random as npr
import numba
from scipy import ndimage

import itertools


def find_local_max(img, d_rad, threshold=1e-15):
    """
    This is effectively a replacement for pkfnd in the matlab/IDL code.

    The output of this function is meant to be feed into :py:func:`~subpixel_centroid`

    The magic of numpy means this should work for any dimension data.

    :param img: an ndarray representing the data to find the local maxes
    :param d_rad: the radius of the dilation, the smallest possible spacing between local maximum
    :param threshold: optional, voxels < threshold are ignored.

    :rtype: (d,N) array of the local maximums.
    """
    d_rad = int(d_rad)
    img = np.array(np.squeeze(img))       # knock out singleton dimensions
    img[img < threshold] = -np.inf        # mask out pixels below threshold
    dim = img.ndim                        # get the dimension of data

    # make structuring element
    s = ndimage.generate_binary_structure(dim, 1)
    # scale it up to the desired size
    d_struct = ndimage.iterate_structure(s, int(d_rad))
    dilated_img = ndimage.grey_dilation(img,
                                        footprint=d_struct,
                                        cval=0,
                                        mode='constant')   # do the dilation

    # find the locations that are the local maximum
    # TODO clean this up
    local_max = np.where(np.exp(img - dilated_img) > (1 - 1e-15))
    # the extra [::-1] is because matplotlib and ndimage disagree an xy vs yx
    return np.vstack(local_max[::-1])

@numba.autojit
def _refine_centroids_loop(img, local_maxes, mask_rad, offset_masks, d_struct, r2_mask):
    results = np.zeros((4, local_maxes.shape[1]), dtype=np.float32)
    for i in range(local_maxes.shape[1]):
        x = local_maxes[1, i]
        y = local_maxes[0, i]
        mass = 0.
        shiftx_accum = 0.
        shifty_accum = 0.
        r2 = 0.
        imd = 0.
        for xi in range(2 * mask_rad + 1):
            for yi in range(2 * mask_rad + 1):
                if d_struct[xi, yi]:
                    imd = img[y + yi - mask_rad, x + xi - mask_rad]
                    mass += imd * d_struct[xi, yi]
                    shiftx_accum += imd * offset_masks[0, xi, yi]
                    shifty_accum += imd * offset_masks[1, xi, yi]
                    r2 += imd * r2_mask[xi, yi]
        results[0, i] = shifty_accum / mass # Note that local_maxes has xy backwards.
        results[1, i] = shiftx_accum / mass
        results[2, i] = mass
        results[3, i] = r2
    return results
def subpixel_centroid(img, local_maxes, mask_rad, struct_shape='circle'):
    '''
    This is effectively a replacement for cntrd in the matlab/IDL code.

    Works for 2D data only. Accelerated by numba.

    :param img: the data
    :param local_maxes: a (d,N) array with the location of the local maximums (as generated by :py:func:`~find_local_max`)
    :param mask_rad: the radius of the mask used for the averaging.
    :param struct_shape: ['circle' | 'diamond'] Shape of mask over each particle.

    :rtype: (d,N) array of positions, (d,) array of masses, (d,) array of r2,
    '''
    local_maxes = local_maxes[::-1]
    # do some data checking/munging
    img = np.squeeze(img)                 # knock out singleton dimensions
    dim = img.ndim
    so = [slice(-mask_rad, mask_rad + 1)] * dim
    # Make circular structuring element
    if struct_shape == 'circle':
        d_struct = (np.sum(np.mgrid[so]**2, 0) < mask_rad**2).astype(np.int8) # NCK IMPROVEMENT
    elif struct_shape == 'diamond':
        s = ndimage.generate_binary_structure(dim, 1)
        # scale it up to the desired size
        d_struct = ndimage.iterate_structure(s, int(mask_rad))
    else: raise ValueError('Shape must be diamond or circle')
    
    offset_masks = np.array([d_struct * os for os in np.mgrid[so]]).astype(np.int8)
    
    r2_mask = np.zeros(d_struct.shape)
    for o in offset_masks:
        r2_mask += o ** 2
    r2_mask = np.sqrt(r2_mask).astype(float)
    results = _refine_centroids_loop(img, local_maxes, mask_rad, offset_masks, d_struct, r2_mask)
    pos = (results[0:2,:] + local_maxes)[::-1,:]
    m = results[2,:]
    r2 = results[3,:]
    return pos, m, r2

def subpixel_centroid_nd(img, local_maxes, mask_rad):
    '''
    This is effectively a replacement for cntrd in the matlab/IDL code.

    Should work for any dimension data


    :param img: the data
    :param local_maxes: a (d,N) array with the location of the local maximums (as generated by :py:func:`~find_local_max`)
    :param mask_rad: the radius of the mask used for the averaging.

    :rtype: (d,N) array of positions, (d,) array of masses, (d,) array of r2,
    '''
    local_maxes = local_maxes[::-1]
    # do some data checking/munging
    mask_rad = int(mask_rad)
    img = np.squeeze(img)                 # knock out singleton dimensions
    # make sure local_maxes.shape makes sense
    dim = img.ndim
    s = ndimage.generate_binary_structure(dim, 1)
    # scale it up to the desired size
    d_struct = ndimage.iterate_structure(s, int(mask_rad))

    so = [slice(-mask_rad, mask_rad + 1)] * dim
    offset_masks = [d_struct * os for os in np.mgrid[so]]

    r2_mask = np.zeros(d_struct.shape)
    for o in offset_masks:
        r2_mask += o ** 2

    r2_mask = np.sqrt(r2_mask)

    shifts_lst = []
    mass_lst = []
    r2_lst = []
    for loc in itertools.izip(*local_maxes):

        window = [slice(p - mask_rad, p + mask_rad + 1) for p in loc]
        img_win = img[window]
        mass = np.sum(img_win * d_struct)
        mass_lst.append(mass)
        shifts_lst.append([np.sum(img_win * o) / mass for o in offset_masks])
        r2_lst.append(np.sum(r2_mask * img_win))
    sub_pixel = np.array(shifts_lst).T + local_maxes
    return sub_pixel[::-1], mass_lst, r2_lst


def band_pass(img, p_rad, hwhm):
    '''
    Intended to be a replacement for bpass in the matlab/IDL code.

    Works by convolving a Gaussian with the image, than a box car and
    taking the difference.

    :param img: array of data
    :param p_rad: the size of the window used for the convolution
    :param hwhm: the hwhm of the Gaussian
    :rtype: :class:`numpy.ndarray` scaled between 0 and 1
    '''
    # make sure the input data is an array and float type.
    img = np.asarray(img).astype(float)

    p_dia = 2 * p_rad + 1

    # do the two convolutions.
    # These should maybe be replaced with masked kernels, but this is
    # faster to code up.
    img_boxcar = ndimage.filters.uniform_filter(img, p_dia, mode='nearest', cval=0)
    img_gaus = ndimage.filters.gaussian_filter(img, hwhm, mode='nearest', cval=0)

    # subtract them
    ret_img = img_boxcar - img_gaus

    # kill data at edegs where the convolution leaked out
    ret_img[ret_img < 0] = 0
    ret_img[:p_dia, :] = 0
    ret_img[-p_dia:, :] = 0
    ret_img[:, :p_dia] = 0
    ret_img[:, -p_dia:] = 0

    # normalize the image
    ret_img -= np.min(ret_img)
    ret_img /= np.max(ret_img)

    return ret_img


def gen_fake_data(list_of_locs, p_rad, hwhm, img_shape):
    """
    Function to generate fake images for testing purposes
    """
    img = np.zeros(img_shape)

    def pixel_values(window, loc):
        i = np.mgrid[window] - loc.reshape(len(window), *[1] * len(window))
        r = np.zeros(i[0].shape)
        for _ in i:
            r += _ ** 2

        return np.exp(-r / (hwhm ** 2))

    for loc in itertools.izip(*list_of_locs):
        window = [slice(int(p) - (p_rad + 2), int(p) + (p_rad + 2) + 1) for p in loc]

        p = pixel_values(window, np.array(loc))
        img[window] += p

    img *= 5
    img += npr.randn(*img.shape) * .1

    return img
