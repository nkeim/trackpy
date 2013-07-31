from .storage import *

import tempfile
import numpy

def test_Pandas():
    data = pandas.Series(numpy.random.random((100,)))
    tmpf = tempfile.NamedTemporaryFile(suffix='.h5', delete=False)
    tmpf.close()
    pobj = Pandas(tmpf.name)
    try:
        pobj.save(data)
        databounced = pobj.read()
        assert numpy.allclose(data.values, databounced.values)
    finally:
        os.unlink(pobj.filename)

def test_JSON():
    data = dict(one=1, two=2.0, three=[4, 5])
    tmpf = tempfile.NamedTemporaryFile(suffix='.json', delete=False)
    tmpf.close()
    pobj = JSON(tmpf.name)
    try:
        pobj.save(data)
        databounced = pobj.read()
        for k in data:
            d = data[k]
            bd = databounced[k]
            assert d == bd
            assert type(d) == type(bd)
    finally:
        os.unlink(pobj.filename)

