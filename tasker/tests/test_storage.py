import os
from path import path

import tempfile, shutil
import unittest
import numpy
import pandas

from tasker.storage import Pandas, Pickle, JSON

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

def test_Pickle():
    data = dict(one=1, two=2.0, three=[4, 5], s={1, 2})
    tmpf = tempfile.NamedTemporaryFile(suffix='.pickle', delete=False)
    tmpf.close()
    pobj = Pickle(tmpf.name)
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

class parentdir(unittest.TestCase):
    def test_parentdir(self):
        """Check whether the file reference can be made absolute after construction."""
        data = dict(one=1, two=2.0, three=[4, 5])
        testdir = path(tempfile.mkdtemp(dir='.'))
        try:
            pobj = JSON('test.json')
            pobj.set_parentdir(testdir)
            pobj.save(data)
            pobj.read()
            assert (testdir / 'test.json').exists()

            # Check auto subdirectory creation (FileBase._mkdir())
            pobj = JSON('subdir/test.json')
            pobj.set_parentdir(testdir)
            pobj.save(data)
            pobj.read()
            assert (testdir / 'subdir' / 'test.json').exists()

            pobj = Pandas('subdir/test.h5')
            pobj.set_parentdir(testdir)
            pobj.save(pandas.DataFrame([data,]))
            pobj.read()
            assert (testdir / 'subdir' / 'test.h5').exists()
        finally:
            shutil.rmtree(testdir)
