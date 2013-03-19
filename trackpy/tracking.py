#Copyright 2012 Thomas A Caswell
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
from scipy.spatial import cKDTree
import numba


class TreeFinder(object):
    def __init__(self, points):
        """Takes a list of particles.
        """
        self.points = points
        coords = np.array([pt.pos for pt in points])
        self.kdtree = cKDTree(coords, int(round(np.log10(len(points)))))
    def get_region(self, point, rrange):
        '''
        :param point: point to find the features around
        :param rrange: the size of the ball to search in


        Returns all the particles with in the region of maximum radius
        rrange in data units
        '''
        # We hard-code the max number of points found.
        # If get_region() returns more particles than this,
        # the search distance is clearly too large.
        dists, inds = self.kdtree.query(point.pos, 10, distance_upper_bound=rrange)
        finite = ~np.isinf(dists)
        return [self.points[i] for i in inds.compress(finite)], dists.compress(finite)

class TreeFinderBare(object):
    def __init__(self, coords):
        """Takes a list of particles.
        """
        self.kdtree = cKDTree(coords, int(round(np.log10(len(coords)))))
class Track(object):
    '''
    :param point: The first feature in the track if not  `None`.
    :type point: :py:class:`~trackpy.tracking.Point`

    Base class for objects to represent linked tracks.  Includes logic
    for adding, removing features to the track.  This can be sub-classed
    to provide additional track level computation as needed.


    '''
    count = 0

    def __init__(self, point=None):
        self.points = []
        # will take initiator point
        if not point is None:
            self.add_point(point)

        self.indx = Track.count           # unique id
        Track.count += 1

    def __iter__(self):
        return self.points.__iter__()

    def __len__(self):
        return len(self.points)

    def __eq__(self, other):
        return self.index == other.index

    def __neq__(self, other):
        return not self.__eq__(other)
    __hash__ = None

    def add_point(self, point):
        '''
        :param point: point to add
        :type point:  :py:class:`~trackpy.tracking.Point`

        Appends the point to this track. '''
        self.points.append(point)
        point.add_to_track(self)

    def remove_point(self, point):
        '''
        :param point: point to remove from this track
        :type point:  :py:class:`~trackpy.tracking.Point`

        removes a point from this track'''
        self.points.remove(point)
        point.track = None

    def last_point(self):
        '''
        :rtype: :py:class:`~trackpy.tracking.Point`

        Returns the last point on the track'''
        return self.points[-1]


class Point(object):
    '''
    Base class for point (features) used in tracking.  This class
    contains all of the general stuff for interacting with
    :py:class:`~trackpy.tracking.Track` objects.



    .. note:: To be used for tracking this class must be sub-classed to provide a :py:func:`distance` function.  Child classes **MUST** call :py:func:`Point.__init__`.  (See :py:class:`~trackpy.tracking.PointND` for example. )
    '''
    count = 0

    def __init__(self):
        self.track = None
        self.uuid = Point.count         # unique id for __hash__
        Point.count += 1

    def __hash__(self):
        return self.uuid

    def __eq__(self, other):
        return self.uuid == other.uuid

    def __neq__(self, other):
        return not self.__eq__(other)

    def add_to_track(self, track):
        '''
        :param track: the track to assign to this :py:class:`Point`

        Sets the track of a :py:class:`Point` object.  Raises
        :py:exc:`Exception` if the object is already assigned a track.



        '''
        if self.track is not None:
            raise Exception("trying to add a particle already in a track")
        self.track = track

    def remove_from_track(self, track):
        '''
        :param track: the track to disassociate from this :py:class:`Point`

        Removes this point from the given track. Raises :py:exc:`Exception` if
        particle not associated with the given track.


        '''
        if self.track != track:
            raise Exception("Point not associated with given track")
        track.remove_point(self)

    def in_track(self):
        '''
        :rtype: bool

        Returns if a point is associated with a track '''
        return self.track is not None


class PointND(Point):
    '''
    :param t: a time-like variable.
    :param pos: position of feature
    :type pos: iterable of length d

    Version of :py:class:`Point` for tracking in flat space with
    non-periodic boundary conditions.
    '''

    def __init__(self, t, pos):
        Point.__init__(self)                  # initialize base class
        self.t = t                            # time
        self.pos = np.asarray(pos)            # position in ND space

    def distance(self, other_point):
        '''
        :param other_point: point to get distance to.
        :type other_point: :py:class:`~trackpy.tracking.Point`

        Returns the absolute distance between this point and other_point

        '''
        return np.sqrt(np.sum((self.pos - other_point.pos) ** 2))

class Level(object):
    def __init__(self, coordinates, levelnum, max_candidates=9):
        self.coordinates = coordinates
        self.levelnum = levelnum
        self.max_candidates = max_candidates
        self.n = len(self.coordinates)
        self.trackassign = np.ones((self.n,), dtype=np.int64) * -1
        self.forward_cands = np.ones((self.n, self.max_candidates), dtype=np.int64) * -1
        self.forward_dists = np.ones((self.n, self.max_candidates), dtype=np.float64) * np.inf
        self.forward_ncands = np.zeros((self.n,), dtype=np.int64)
        self.backward_cands = np.ones_like(self.forward_cands) * -1
        self.backward_dists = np.ones_like(self.forward_dists) * np.inf
        self.backward_ncands = np.zeros_like(self.forward_ncands)
    def add_mem_points(self):
        """Add extra points from the memory buffer"""
        raise NotImplementedError()
def fastlink(levels, search_range, memory=0, track_cls=Track, diagdict={}):
    '''
    :param levels: Iterable of (levelnum, coords) tuples
    :param search_range: the maximum distance features can move between frames
    :param memory: the maximum number of frames that a feature can skip along a track
    :param track_cls: The class to use for the returned track objects
    :type track_cls: :py:class:`~trackpy.tracking.Track`

    Generic version of linking algorithm, should work for any
    dimension.  All information about dimensionality and the metric
    are encapsulated in the hash_table and
    :py:class:`~trackpy.tracking.Point` objects.

    '''
    #    print "starting linking"
    # initial source set
    lev_iter = iter(levels)
    prev_levelnum, prev_levelcoords = lev_iter.next()
    prev_level = Level(prev_levelcoords, prev_levelnum)
    prev_set = set(range(prev_level.n))
    prev_hash = TreeFinderBare(prev_level.coordinates)

    # assume everything in first level starts a track
    prev_level.trackassign = np.arange(prev_level.n, dtype=prev_level.trackassign.dtype)
    yield prev_level.trackassign
    maxtrack = prev_level.n - 1
    #mem_set = set()
    # fill in first 'prev' hash

    # fill in memory list of sets
    #mem_history = []
    #for j in range(memory):
        #mem_history.append(set())

    for cur_levelnum, cur_levelcoords in lev_iter:
        cur_level = Level(cur_levelcoords, cur_levelnum)
        # create the set for the destination level
        cur_set = set(range(cur_level.n))
        # make a new hash object
        cur_hash = TreeFinderBare(cur_level.coordinates)

        # create a second copy that will be used as the source in
        # the next loop
        tmp_set = set(range(cur_level.n))
        # memory set
        #new_mem_set = set()

        # Look for particles in the previous level near a current-level particle's position
        # These are our backward candidates
        maxcands = cur_level.max_candidates
        for i in xrange(len(cur_level.coordinates)):
            dists, inds = prev_hash.kdtree.query(cur_level.coordinates[i], maxcands, distance_upper_bound=search_range)
            cur_level.backward_cands[i] = inds
            cur_level.backward_dists[i] = dists
        # Invert cur_level.backward into prev_level.forward
        _invert_candidates(cur_level.backward_cands, cur_level.backward_dists,
                prev_level.forward_cands, prev_level.forward_dists)
        # Count up the number of candidates for each particle
        cur_level.backward_ncands = np.sum(~np.isinf(cur_level.backward_dists), axis=1)
        prev_level.forward_ncands = np.sum(~np.isinf(prev_level.forward_dists), axis=1)

        # Do the trivial bonds all at once
        # First, ones with no backwards candidates. They each need a new track ID
        singles = np.where(cur_level.backward_ncands == 0)[0]
        cur_level.trackassign.put(singles, np.arange(maxtrack+1, maxtrack+len(singles)+1,
            dtype=cur_level.trackassign.dtype))
        maxtrack += len(singles)
        # Remove these particles from consideration
        cur_set.difference_update(set(singles.tolist()))

        # Next, where the only backwards candidate also has only one forwards candidate
        # This has to be done in multiple steps. numpy needs a join() operation like pandas!
        couple_cands_cpid = np.where(cur_level.backward_ncands == 1)[0]
        # The ID of that one old candidate should be in the first column
        couple_cands_ppid = cur_level.backward_cands[couple_cands_cpid,0] 
        # Find out how many candidates those old particles have
        couple_cands_fwd_iscouple = prev_level.forward_ncands.take(couple_cands_ppid) == 1
        # Winnow down the IDs to just those couples
        couples_cpid = couple_cands_cpid.compress(couple_cands_fwd_iscouple)
        couples_ppid = couple_cands_ppid.compress(couple_cands_fwd_iscouple)
        # Copy the track assignments to the new level
        cur_level.trackassign.put(couples_cpid, prev_level.trackassign.take(couples_ppid))
        # Remove these particles from consideration
        cur_set.difference_update(set(couple_cands_cpid.tolist()))
        prev_set.difference_update(set(couple_cands_ppid.tolist()))

        #diagdict['sfcands'] = prev_level.forward_cands
        #diagdict['sfdists'] = prev_level.forward_dists
        #diagdict['dbcands'] = cur_level.backward_cands
        #diagdict['dbdists'] = cur_level.backward_dists

        # Now we are ready to do everyone else!
        # while there are particles left to link, link
        while len(prev_set) > 0 and len(cur_set) > 0:
            cpid = cur_set.pop()
            #bc_c = cur_level.backward_ncands[cpid]
            # we need to generate the sub networks
            done_flg = False
            s_sn = set()                  # source sub net
            d_sn = set()                  # destination sub net
            # add working particle to destination sub-net
            d_sn.add(cpid)
            while not done_flg:
                d_sn_sz = len(d_sn)
                s_sn_sz = len(s_sn)
                for dp in d_sn:
                    for c_sp in cur_level.backward_cands[dp,:cur_level.backward_ncands[dp]]:
                        s_sn.add(c_sp)
                        prev_set.discard(c_sp)
                for sp in s_sn:
                    for c_dp in prev_level.forward_cands[sp,:prev_level.forward_ncands[sp]]:
                        d_sn.add(c_dp)
                        cur_set.discard(c_dp)
                done_flg = (len(d_sn) == d_sn_sz) and (len(s_sn) == s_sn_sz)

            s_sn_list = sorted(s_sn)
            best_pairs = link_subnet_bare(s_sn_list, prev_level.forward_cands[s_sn_list,:],
                    prev_level.forward_dists[s_sn_list,:],
                    prev_level.forward_ncands[s_sn_list],
                    search_range)
            spl, dpl = zip(*best_pairs)

            # remove the linked particles
            d_remain = d_sn - set(dpl)
            #s_remain = s_sn - set(spl)

            # do linking
            cur_level.trackassign.put(dpl, prev_level.trackassign.take(spl))
            for dp in d_remain:
                # if unclaimed destination particle, a track is born!
                maxtrack += 1
                cur_level.trackassign[dp] = maxtrack
            # tack all of the unmatched source particles into the new
            # memory set
            #new_mem_set |= s_remain

        # set prev_hash to cur hash
        prev_hash = cur_hash
        # add in the memory points
        # store the current level for use in next loop
#        if memory > 0:
#            # identify the new memory points
#            new_mem_set -= mem_set
#            mem_history.append(new_mem_set)
#            # remove points that are now too old
#            mem_set -= mem_history.pop(0)
#            # add the new points
#            mem_set |= new_mem_set
#            # add the memory particles to what will be the next source
#            # set
#            tmp_set |= mem_set # add memory points to prev_hash (to be used as the next source)
#            for m in mem_set:
#                # add points to the hash
#                prev_hash.add_point(m)
#                # re-create the forward_cands list
#                m.forward_cands = []
        prev_set = tmp_set
        prev_level = cur_level

        # Some particles were denied tracks by the subnet code. They need new ones!
        loners = (cur_level.trackassign < 0).nonzero()[0]
        cur_level.trackassign.put(loners, np.arange(maxtrack+1, maxtrack+len(loners)+1,
            dtype=cur_level.trackassign.dtype))
        maxtrack += len(loners)

        # No particles unaccounted for
        assert len((cur_level.trackassign < 0).nonzero()[0]) == 0 
        yield cur_level.trackassign


def link(levels, search_range, memory=0, track_cls=Track):
    '''
    :param levels: Nested iterables of :py:class:`~trapy.tracking.Point` objects
    :type levels: Iterable of iterables
    :param search_range: the maximum distance features can move between frames
    :param memory: the maximum number of frames that a feature can skip along a track
    :param track_cls: The class to use for the returned track objects
    :type track_cls: :py:class:`~trackpy.tracking.Track`

    Generic version of linking algorithm, should work for any
    dimension.  All information about dimensionality and the metric
    are encapsulated in the hash_table and
    :py:class:`~trackpy.tracking.Point` objects.

    '''
    #    print "starting linking"
    # initial source set
    lev_iter = iter(levels)
    prev_level = lev_iter.next()
    prev_set = set(prev_level)
    prev_hash = TreeFinder(prev_level)
    # set up the particles in the previous level for
    # linking
    for p in prev_set:
        p.forward_cands = []

    # assume everything in first level starts a track
    # initialize the master track list with the points in the first level
    track_lst = [track_cls(p) for p in prev_set]
    mem_set = set()
    # fill in first 'prev' hash

    # fill in memory list of sets
    mem_history = []
    for j in range(memory):
        mem_history.append(set())

    for cur_level in lev_iter:
        # make a new hash object
        cur_hash = TreeFinder(cur_level)

        # create the set for the destination level
        cur_set = set(cur_level)
        # create a second copy that will be used as the source in
        # the next loop
        tmp_set = set(cur_level)
        # memory set
        new_mem_set = set()

        # fill in first 'cur' hash and set up attributes for keeping
        # track of possible connections

        for p in cur_set:
            p.back_cands = []
            p.forward_cands = []
        # sort out what can go to what
        for p in cur_level:
            # get
            candidates, distances = prev_hash.get_region(p, search_range)
            for wp, d in zip(candidates, distances):
                # FIXME: PointND.distance() is now unused
                p.back_cands.append((wp, d))
                wp.forward_cands.append((p, d))

        # sort the candidate lists by distance
        for p in cur_set:
            p.back_cands.sort(key=lambda x: x[1])
        for p in prev_set:
            p.forward_cands.sort(key=lambda x: x[1])
        # while there are particles left to link, link
        while len(prev_set) > 0 and len(cur_set) > 0:
            p = cur_set.pop()
            bc_c = len(p.back_cands)
            # no backwards candidates
            if bc_c == 0:
                # add a new track
                track_lst.append(track_cls(p))
                # clean up tracking apparatus
                del p.back_cands
                # short circuit loop
                continue
            if bc_c == 1:
                # one backwards candidate
                b_c_p = p.back_cands[0]
                # and only one forward candidate
                if len(b_c_p[0].forward_cands) == 1:
                    # add to the track of the candidate
                    b_c_p[0].track.add_point(p)
                    # clean up tracking apparatus
                    del p.back_cands
                    del b_c_p[0].forward_cands
                    # short circuit loop
                    continue
            # we need to generate the sub networks
            done_flg = False
            s_sn = set()                  # source sub net
            d_sn = set()                  # destination sub net
            # add working particle to destination sub-net
            d_sn.add(p)
            while not done_flg:
                d_sn_sz = len(d_sn)
                s_sn_sz = len(s_sn)
                for dp in d_sn:
                    for c_sp in dp.back_cands:
                        s_sn.add(c_sp[0])
                        prev_set.discard(c_sp[0])
                for sp in s_sn:
                    for c_dp in sp.forward_cands:
                        d_sn.add(c_dp[0])
                        cur_set.discard(c_dp[0])
                done_flg = (len(d_sn) == d_sn_sz) and (len(s_sn) == s_sn_sz)

            #snl = sub_net_linker(s_sn, search_range)

            #spl, dpl = zip(*snl.best_pairs)
            best_pairs = link_subnet(s_sn, search_range)
            spl, dpl = zip(*best_pairs)


            # strip the distance information off the subnet sets and
            # remove the linked particles
            d_remain = set([d for d in d_sn])
            d_remain -= set(dpl)
            s_remain = set([s for s in s_sn])
            s_remain -= set(spl)

            for sp, dp in best_pairs:
                # do linking and clean up
                sp.track.add_point(dp)
                del dp.back_cands
                del sp.forward_cands
            for sp in s_remain:
                # clean up
                del sp.forward_cands
            for dp in d_remain:
                # if unclaimed destination particle, a track in born!
                track_lst.append(track_cls(dp))
                # clean up
                del dp.back_cands
            # tack all of the unmatched source particles into the new
            # memory set
            new_mem_set |= s_remain

        # set prev_hash to cur hash
        prev_hash = cur_hash
        if memory > 0:
            # identify the new memory points
            new_mem_set -= mem_set
            mem_history.append(new_mem_set)
            # remove points that are now too old
            mem_set -= mem_history.pop(0)
            # add the new points
            mem_set |= new_mem_set
            # add the memory particles to what will be the next source
            # set
            tmp_set |= mem_set # add memory points to prev_hash (to be used as the next source)
            for m in mem_set:
                # add points to the hash
                prev_hash.add_point(m)
                # re-create the forward_cands list
                m.forward_cands = []
        prev_set = tmp_set

        # add in the memory points
        # store the current level for use in next loop

    return track_lst

@numba.autojit
def _invert_candidates(back_cands, back_dists, fwd_cands, fwd_dists):
    """Take the information in back_* and reshuffle into fwd_*.

    Modifies fwd_*; returns nothing.
    """
    inf = np.inf
    maxcands = fwd_cands.shape[1]
    for i in range(back_cands.shape[0]):
        for j in range(maxcands):
            ppid = back_cands[i,j]
            pdist = back_dists[i,j]
            if pdist != inf:
                for fj in range(maxcands):
                    if fwd_dists[ppid,fj] == inf:
                        # Find an empty spot in fwd_cands
                        fwd_dists[ppid,fj] = pdist
                        fwd_cands[ppid,fj] = i
                        break
    return 0

#######
# Sub-network optimization code

class SubnetOversizeException(Exception):
    '''An :py:exc:`Exception` to be raised when the sub-nets are too
    big to be efficiently linked.  If you get this then either reduce your search range
    or increase :py:attr:`sub_net_linker.MAX_SUB_NET_SIZE`'''
    pass

def link_subnet(s_sn, search_radius):
    """Recursively find the optimal bonds for a group of particles between 2 frames.
    
    This is only invoked when there is more than one possibility within
    'search_radius'.
    """
    # The basic idea: replace Point objects with integer indices into lists of Points.
    # Then the hard part (recursion) runs quickly because it is just passing arrays.
    # In fact, we can compile it with numba so that it runs in acceptable time.
    MAX_SUB_NET_SIZE = 50 # Can't exceed Python's recursion depth
    max_candidates = 9 # Max forward candidates we expect for any particle
    src_net = list(s_sn)
    nj = len(src_net) # j will index the source particles
    if nj > MAX_SUB_NET_SIZE:
        raise SubnetOversizeException('sub net contains %d points' % nj)
    # Build arrays of all destination (forward) candidates and their distances
    dcands = set()
    for p in src_net:
        dcands.update([cand for cand, dist in p.forward_cands])
    dcands = list(dcands)
    dcands_map = {cand: i for i, cand in enumerate(dcands)}
    # A source particle's actual candidates only take up the start of
    # each row of the array. All other elements represent the null link option
    # (i.e. particle lost)
    candsarray = np.ones((nj, max_candidates + 1), dtype=np.int64) * -1
    distsarray = np.ones((nj, max_candidates + 1), dtype=np.float64) * search_radius
    ncands = np.zeros((nj,), dtype=np.int64)
    for j, sp in enumerate(src_net):
        ncands[j] = len(sp.forward_cands)
        if ncands[j] > max_candidates:
            raise SubnetOversizeException('Particle has %i forward candidates --- too many' % ncands[j])
        candsarray[j,:ncands[j]] = [dcands_map[cand] for cand, dist in sp.forward_cands]
        distsarray[j,:ncands[j]] = [dist for cand, dist in sp.forward_cands]
    # The assignments are persistent across levels of the recursion
    best_assignments = np.ones((nj,), dtype=np.int64) * -1
    cur_assignments = np.ones((nj,), dtype=np.int64) * -1
    snr = _SNRecursion()
    snr.sn_recur(0, nj, np.inf, 0., ncands, candsarray, distsarray, best_assignments, cur_assignments)
    # Remove null links and return particle objects
    return [(src_net[j], dcands[i]) for j, i in enumerate(best_assignments) if i >= 0]

def link_subnet_bare(subnet_src_pids, subnet_fwd_cands, subnet_fwd_dists, 
        subnet_fwd_ncands, search_radius):
    """Recursively find the optimal bonds for a group of particles between 2 frames.
    
    This is only invoked when there is more than one possibility within
    'search_radius'.
    """
    # The basic idea: replace Point objects with integer indices into lists of Points.
    # Then the hard part (recursion) runs quickly because it is just passing arrays.
    # In fact, we can compile it with numba so that it runs in acceptable time.
    #print 'subnet'
    MAX_SUB_NET_SIZE = 50 # Can't exceed Python's recursion depth
    src_net = list(subnet_src_pids)
    nj = len(src_net) # j will index the source particles
    if nj > MAX_SUB_NET_SIZE:
        raise SubnetOversizeException('sub net contains %d points' % nj)
    # Build arrays of all destination (forward) candidates and their distances
    # A source particle's actual candidates only take up the start of
    # each row of the array. All other elements represent the null link option
    # (i.e. particle lost)
    # Conveniently, this is just like the format of the candidate arrays
    # we were passed from fastlink(), except that the null option has
    # distance 'search_radius', not infinity.
    subnet_fwd_dists_sr = np.clip(subnet_fwd_dists, 0, search_radius)
    # We also need to sort by distance, so we try the most likely assignments first
    sorter = np.argsort(subnet_fwd_dists_sr, axis=1)
    sorthelper = np.mgrid[0:sorter.shape[0],0:sorter.shape[1]][0]
    sfdists = subnet_fwd_dists_sr[sorthelper,sorter]
    sfcands = subnet_fwd_cands[sorthelper,sorter]
    # FIXME: Technically, we need an extra null-link column, in case a particle
    # has the maximum number of candidates.
    # But in that case, the subnet code will take so long anyhow that the user will
    # try a lower search_radius :)

    # The assignments are persistent across levels of the recursion
    best_assignments = np.ones((nj,), dtype=np.int64) * -1
    cur_assignments = np.ones((nj,), dtype=np.int64) * -1
    tmp_assignments = np.zeros((nj,), dtype=np.int64)
    # FIXME: This can be done without recursion! We just have to implement
    # our own "stack" with 'nj' levels. That would be a massive speedup.
    #print  subnet_fwd_ncands, sfcands, sfdists
    #snr = _SNRecursion()
    #snr.sn_recur(0, nj, np.inf, 0., subnet_fwd_ncands, sfcands, 
    #        sfdists, best_assignments, cur_assignments)
    _sn_norecur(subnet_fwd_ncands, sfcands, sfdists, cur_assignments, tmp_assignments, best_assignments)
    #print best_assignments
    # Remove null links and return particle objects
    return [(src_net[j], i) for j, i in enumerate(best_assignments) if i >= 0]

@numba.jit
class _SNRecursion(object):
    @numba.f8(numba.i8, numba.i8, numba.f8, numba.f8,
        numba.i8[:],
        numba.i8[:,:], numba.f8[:,:], numba.i8[:], numba.i8[:])
    def sn_recur(self, j, nj, best_sum, cur_sum, ncands, candsarray, distsarray, best_assignments, cur_assignments):
        # ncands[j] is number of non-nan elements in candsarray[j]
        for i in range(ncands[j] + 1): # Include the null link
            tmp_sum = cur_sum + distsarray[j,i]
            if tmp_sum > best_sum:
                # if we are already greater than the best sum, bail. we
                # can bail all the way out of this branch because all
                # the other possible connections (including the null
                # connection) are more expensive than the current
                # connection, thus we can discard with out testing all
                # leaves down this branch
                return best_sum
            # We can have as many null links as we want, but the real particles are finite
            flag = 0
            for jtmp in range(nj): 
                if cur_assignments[jtmp] == candsarray[j,i]:
                    flag = 1
            if flag and candsarray[j,i] >= 0:
                # we have already used this destination point, bail
                continue
            # OK, I guess we'll try this assignment
            cur_assignments[j] = candsarray[j,i]
            if j + 1 == nj:
                # We have made assignments for all the particles,
                # and we never exceeded the previous best_sum.
                # This is our new optimum.
                #print 'hit: %f' % best_sum
                best_sum = tmp_sum
                # This array is shared by all levels of recursion.
                # If it's not touched again, it will be used once we
                # get back to link_subnet
                for tmpj in range(nj):
                    best_assignments[tmpj] = cur_assignments[tmpj]
            else:
                # Try various assignments for the next particle
                best_sum = self.sn_recur(j + 1, nj, best_sum, tmp_sum, ncands, candsarray, distsarray, best_assignments, cur_assignments)
            # Undo the assignment, for sanity.
            cur_assignments[j] = -1
        return best_sum

@numba.autojit
def _sn_norecur(ncands, candsarray, distsarray, cur_assignments, tmp_assignments, best_assignments):
    """Find the optimal track assigments for a subnetwork, without recursion.

    This is for nj source particles. All arguments are arrays with nj rows.

    cur_assignments, tmp_assignments are just temporary registers of length 
    best_assignments is modified in place.
    Returns the best sum.
    """
    nj = candsarray.shape[0]
    tmp_sum = 0.
    cur_sum = 0.
    best_sum = 1.0e23
    j = 0
    while 1:
        # This is an endless loop. We go up and down levels of recursion,
        # and emulate the mechanics of nested "for" loops, using the
        # blocks of code marked "GO UP" and "GO DOWN". It's not pretty.

        # Load state from the "stack"
        i = tmp_assignments[j]
        if i > ncands[j]:
            # We've exhausted possibilities at this level, including the
            # null link; make no more changes and go up a level
            #### GO UP
            if j > 0:
                j -= 1 
                tmp_assignments[j] += 1 # Try the next candidate at this higher level
                continue
            else:
                return best_sum
        tmp_sum = cur_sum + distsarray[j,i]
        if tmp_sum > best_sum:
            # if we are already greater than the best sum, bail. we
            # can bail all the way out of this branch because all
            # the other possible connections (including the null
            # connection) are more expensive than the current
            # connection, thus we can discard with out testing all
            # leaves down this branch
            #### GO UP
            if j > 0:
                j -= 1
                tmp_assignments[j] += 1 # Try the next candidate at this higher level
                continue
            else:
                return best_sum
        # We can have as many null links as we want, but the real particles are finite
        flag = 0
        for jtmp in range(nj): 
            if cur_assignments[jtmp] == candsarray[j,i]:
                flag = 1
        if flag and candsarray[j,i] >= 0:
            # we have already used this destination point; try the next one instead
            tmp_assignments[j] += 1
            continue
        # OK, I guess we'll try this assignment
        cur_assignments[j] = candsarray[j,i]
        if j + 1 == nj:
            # We have made assignments for all the particles,
            # and we never exceeded the previous best_sum.
            # This is our new optimum.
            #print 'hit: %f' % best_sum
            best_sum = tmp_sum
            # This array is shared by all levels of recursion.
            # If it's not touched again, it will be used once we
            # get back to link_subnet
            for tmpj in range(nj):
                best_assignments[tmpj] = cur_assignments[tmpj]
            #### GO UP
            if j > 0:
                j -= 1 
                tmp_assignments[j] += 1 # Try the next candidate at this higher level
                continue
            else:
                return best_sum
        else:
            # Try various assignments for the next particle
            #### GO DOWN
            j += 1 
            tmp_assignments[j] = 0





