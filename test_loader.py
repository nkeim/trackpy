import os
import tempfile, unittest
from path import path
from . import loader

basedir = path.getcwd()
mypath = path(__file__)
sample_taskfile = mypath.dirname() / 'sample_taskfile.py'

class loadertests(unittest.TestCase):
    def setUp(self):
        os.chdir(basedir)
        self.testdir = path(tempfile.mkdtemp())
        self.taskfile = self.testdir / 'taskfile.py'
        sample_taskfile.copy(self.taskfile)
        self.temptaskdir = self.testdir / 'top' / 'middle' / 'bottom'
        self.temptaskdir.makedirs_p()
    def tearDown(self):
        os.chdir(basedir)
        self.testdir.rmtree()
    def test_loader(self):
        task = loader.use(self.temptaskdir)
        self.assertEqual(task.two()[0], 2.0)
    def test_loaderfail(self):
        self.taskfile.unlink()
        with self.assertRaises(IOError):
            loader.use(self.temptaskdir)

