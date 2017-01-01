from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import six
import functools
import unittest

import nose.tools
import numpy as np
import pandas

import trackpy
from trackpy import predict

# Catch attempts to set values on an inadvertent copy of a Pandas object.
trackpy.utils.make_pandas_strict()

def mkframe(n=1, Nside=3):
    xg, yg = np.mgrid[:Nside,:Nside]
    dx = (n - 1)
    dy = -(n - 1)
    return pandas.DataFrame(
            dict(x=xg.flatten() + dx, y=yg.flatten() + dy, frame=n))

def link(frames, linker, *args, **kw):
    defaults = {}
    defaults.update(kw)
    return pandas.concat(linker(frames,  *args, **defaults),
                       ignore_index=True)

def get_linked_lengths(frames, linker, *args, **kw):
    """Track particles and return the length of each trajectory."""
    return link(frames, linker, *args, **kw).groupby('particle').x.count()

Nside_oversize = int(np.sqrt(100)) # Make subnet linker fail

class LinkWithPrediction(object):
    def get_linker_iter(self, pred):
        """Return the function that will perform linking with predictor 'pred'.
        This lets us test multiple linking implementations.
        """
        return pred.link_df_iter

    def get_linker(self, pred):
        """Return the function that will perform linking with predictor 'pred'.
        This lets us test multiple linking implementations.
        """
        return pred.link_df

    def get_unwrapped_linker(self):
        return trackpy.link_df_iter


class BaselinePredictTests(LinkWithPrediction, unittest.TestCase):
    def test_null_predict(self):
        """Make sure that a prediction of no motion does not interfere
        with normal tracking.
        """
        # link_df_iter
        pred = predict.NullPredict()
        ll = get_linked_lengths((mkframe(0), mkframe(0.25)),
                                self.get_linker_iter(pred), 0.45)
        assert all(ll.values == 2)

        # link_df
        pred = predict.NullPredict()
        ll_df = self.get_linker(pred)(pandas.concat((mkframe(0),
                                                     mkframe(0.25))), 0.45)
        assert all(ll_df.groupby('particle').x.count().values == 2)

        # Make sure that keyword options are handled correctly
        # (This checks both link_df and link_df_iter)
        features = pandas.concat((mkframe(0), mkframe(0.25)))
        features.rename(columns=lambda n: n + '_', inplace=True)
        pred = predict.NullPredict()
        ll_df = self.get_linker(pred)(features, 0.45, t_column='frame_',
                                      pos_columns=['x_', 'y_'])
        assert all(ll_df.groupby('particle').x_.count().values == 2)

    def test_predict_decorator(self):
        """Make sure that a prediction of no motion does not interfere
        with normal tracking.
        """
        pred = predict.null_predict
        pred_link = functools.partial(self.get_unwrapped_linker(), predictor=pred)
        ll = get_linked_lengths((mkframe(0), mkframe(0.25)),
                                pred_link, 0.45)
        assert all(ll.values == 2)

    def test_fail_predict(self):
        ll = get_linked_lengths((mkframe(0), mkframe(0.25), mkframe(0.65)),
                                self.get_unwrapped_linker(), 0.45)
        assert not all(ll.values == 2)

    @nose.tools.raises(trackpy.SubnetOversizeException)
    def test_subnet_fail(self):
        Nside = Nside_oversize
        ll = get_linked_lengths((mkframe(0, Nside),
                                 mkframe(0.25, Nside),
                                 mkframe(0.75, Nside)), self.get_unwrapped_linker(), 1)

class VelocityPredictTests(LinkWithPrediction):
    def test_simple_predict(self):
        pred = self.predict_class()
        ll = get_linked_lengths((self.mkframe(0), self.mkframe(0.25), self.mkframe(0.65)),
                                self.get_linker_iter(pred), 0.45)
        assert all(ll.values == 3)

    def test_big_predict(self):
        Nside = Nside_oversize
        pred = self.predict_class()
        ll = get_linked_lengths((self.mkframe(0, Nside), self.mkframe(0.25, Nside),
                                 self.mkframe(0.75, Nside)),
                                self.get_linker_iter(pred), 0.45)
        assert all(ll.values == 3)

    def test_predict_newparticles(self):
        """New particles should get the velocities of nearby ones."""
        pred = self.predict_class()
        ll = get_linked_lengths((self.mkframe(0), self.mkframe(0.25), self.mkframe(0.65, 4),
                                self.mkframe(1.05, 4)),
                                self.get_linker_iter(pred), 0.45)
        assert not any(ll.values == 1)

    def test_predict_memory(self):
        pred = self.predict_class()
        frames = [self.mkframe(0), self.mkframe(0.25), self.mkframe(0.65),
                  self.mkframe(1.05), self.mkframe(1.45)]
        ll = get_linked_lengths(frames, self.get_linker_iter(pred), 0.45)
        assert all(ll.values == len(frames))

        # Knock out a particle. Make sure tracking fails.
        frames[3].x[5] = np.nan
        frames[3] = frames[3].dropna()
        pred = self.predict_class()
        tr = link(frames, self.get_linker_iter(pred), 0.45)
        starts = tr.groupby('particle').frame.min()
        ends = tr.groupby('particle').frame.max()
        assert not all(ends - starts == 1.45)

        pred = self.predict_class()
        tr = link(frames, self.get_linker_iter(pred), 0.45, memory=1)
        starts = tr.groupby('particle').frame.min()
        ends = tr.groupby('particle').frame.max()
        assert all(ends - starts == 1.45), 'Prediction with memory fails.'

    def test_predict_diagnostics(self):
        """Minimally test predictor instrumentation."""
        pred = self.instrumented_predict_class()
        Nside = Nside_oversize
        frames = (self.mkframe(0, Nside), self.mkframe(0.25, Nside),
                  self.mkframe(0.75, Nside))
        ll = get_linked_lengths(frames, self.get_linker_iter(pred), 0.45)
        assert all(ll.values == 3)
        diags = pred.dump()
        assert len(diags) == 2
        for i, d in enumerate(diags):
            assert d['t1'] == frames[i+1].frame.iloc[0]
            assert 'state' in d
            assert np.all(d['particledf']['x_act'] == frames[i+1].x)


class NearestVelocityPredictTests(VelocityPredictTests, unittest.TestCase):
    def setUp(self):
        self.predict_class = predict.NearestVelocityPredict
        self.instrumented_predict_class = \
            predict.instrumented()(self.predict_class)
        self.mkframe = mkframe
    def test_initial_guess(self):
        """When an accurate initial velocity is given, velocities
        in the first pair of frames may be large."""
        pred = self.predict_class(
            initial_guess_positions=[(0., 0.)],
            initial_guess_vels=[(1., -1.)])
        ll = get_linked_lengths((self.mkframe(0), self.mkframe(1.),
                                 self.mkframe(2.)),
                                self.get_linker_iter(pred), 0.45)
        assert all(ll.values == 3)

class DriftPredictTests(VelocityPredictTests, unittest.TestCase):
    def setUp(self):
        self.predict_class = predict.DriftPredict
        self.instrumented_predict_class = \
            predict.instrumented()(self.predict_class)
        self.mkframe = mkframe
    def test_initial_guess(self):
        """When an accurate initial velocity is given, velocities
        in the first pair of frames may be large."""
        pred = self.predict_class(initial_guess=(1., -1.))
        ll = get_linked_lengths((self.mkframe(0), self.mkframe(1.),
                                 self.mkframe(2.)),
                                self.get_linker_iter(pred), 0.45)
        assert all(ll.values == 3)


class ChannelPredictXTests(VelocityPredictTests, unittest.TestCase):
    def setUp(self):
        self.predict_class = functools.partial(
            predict.ChannelPredict,  3, minsamples=3)
        self.instrumented_predict_class = functools.partial(
            predict.instrumented()(predict.ChannelPredict), 3, minsamples=3)
        self.mkframe = self._channel_frame
    def _channel_frame(self, n=1, Nside=3):
        xg, yg = np.mgrid[:Nside,:Nside]
        dx = (n - 1) * np.sqrt(2)
        dy = 0.
        return pandas.DataFrame(
                dict(x=xg.flatten() + dx, y=yg.flatten() + dy, frame=n))
    def test_initial_guess(self):
        """When an accurate initial velocity profile is given, velocities
        in the first pair of frames may be large."""
        def _shear_frame(t=1., Nside=4):
            xg, yg = np.mgrid[:Nside,:Nside]
            dx = 0.45 * t * yg
            return pandas.DataFrame(
                dict(x=(xg + dx).flatten(), y=yg.flatten(), frame=t))
        inity = np.arange(4)
        initprof = np.vstack((inity, inity*0.45)).T
        # We need a weird bin size (1.1) to avoid bin boundaries coinciding
        # with particle positions.
        pred = predict.ChannelPredict(1.1, minsamples=3,
                                      initial_profile_guess=initprof)
        print(_shear_frame(1.))
        ll = get_linked_lengths((_shear_frame(0), _shear_frame(1.),
                                 _shear_frame(2), _shear_frame(3)),
                        self.get_linker_iter(pred), 0.45)
        assert all(ll.values == 4)

class ChannelPredictYTests(VelocityPredictTests, unittest.TestCase):
    def setUp(self):
        self.predict_class = functools.partial(
            predict.ChannelPredict, 3, 'y', minsamples=3)
        self.instrumented_predict_class = functools.partial(
            predict.instrumented()(predict.ChannelPredict),
            3, 'y', minsamples=3)
        self.mkframe = self._channel_frame
    def _channel_frame(self, n=1, Nside=3):
        xg, yg = np.mgrid[:Nside,:Nside]
        dx = 0.
        dy = (n - 1) * np.sqrt(2)
        return pandas.DataFrame(
            dict(x=xg.flatten() + dx, y=yg.flatten() + dy, frame=n))


## find_link prediction tests
# Define a mixin that converts a normal prediction test class into one
# that uses find_link.

class FindLinkPredictTest(object):
    def get_linker_iter(self, pred):
        return functools.partial(pred.wrap, trackpy.link_simple_df_iter)

    def get_linker(self, pred):
        return functools.partial(pred.wrap_single, trackpy.link_simple_df_iter)

    def get_unwrapped_linker(self):
        return trackpy.link_simple_df_iter


class FLBaselinePredictTests(FindLinkPredictTest, BaselinePredictTests):
    pass


class FLNearestVelocityPredictTests(FindLinkPredictTest, NearestVelocityPredictTests):
    pass


class FLDriftPredictTests(FindLinkPredictTest, DriftPredictTests):
    pass


class FLChannelPredictXTests(FindLinkPredictTest, ChannelPredictXTests):
    pass


class FLChannelPredictYTests(FindLinkPredictTest, ChannelPredictYTests):
    pass


if __name__ == '__main__':
    unittest.main()
