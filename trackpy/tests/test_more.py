import numpy as np
from trackpy.tracking import PointND, link
from copy import deepcopy

# Call lambda function for a fresh copy each time.
unit_steps = lambda: [[PointND(t, (x, 0))] for t, x in enumerate(range(5))]

np.random.seed(0)
random_x = np.random.randn(5).cumsum()
random_x -= random_x.min()  # All x > 0
max_disp = np.diff(random_x).max()
random_walk = lambda: [[PointND(t, (x, 5))] for t, x in enumerate(random_x)]


def test_search_range():
    t = link(unit_steps(), 1.1)
    assert len(t) == 1  # One track
    t_short = link(unit_steps(), 0.9)
    assert len(t_short) == len(unit_steps())  # Each step is a separate track.

    t = link(random_walk(), max_disp + 0.1)
    assert len(t) == 1  # One track
    t_short = link(random_walk(), max_disp - 0.1)
    assert len(t_short) > 1  # Multiple tracks


def test_memory():
    """A unit-stepping trajectory and a random walk are observed
    simultaneously. The random walk is missing from one observation."""
    a = [p[0] for p in unit_steps()]
    b = [p[0] for p in random_walk()]
    # b[2] is intentionally omitted below.
    gapped = lambda: deepcopy([[a[0], b[0]], [a[1], b[1]], [a[2]],
                              [a[3], b[3]], [a[4], b[4]]])
    safe_disp = 1 + random_x.max() - random_x.min()  # Definitely large enough
    t0 = link(gapped(), safe_disp, memory=0)
    assert len(t0) == 3, len(t0)
    t5 = link(gapped(), safe_disp, memory=5)
    assert len(t5) == 2, len(t5)

