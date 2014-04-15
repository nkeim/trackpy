import os
import tempfile, unittest
from path import path
from tasker import loader

basedir = path.getcwd()
mypath = path(__file__)
sample_taskfile = mypath.dirname() / 'sample_taskfile.py'

class loadertests(unittest.TestCase):
    def setUp(self):
        os.chdir(basedir)
        self.testdir = path(tempfile.mkdtemp())
        self.taskfile = self.testdir / 'taskfile_sub.py'
        sample_taskfile.copy(self.taskfile)
        self.temptaskdir = self.testdir / 'top' / 'middle' / 'bottom'
        self.temptaskdir.makedirs_p()
    def tearDown(self):
        os.chdir(basedir)
        self.testdir.rmtree()
    def test_loader(self):
        assert loader.which(self.temptaskdir).basename() == 'taskfile_sub.py'
        task = loader.use(self.temptaskdir)
        self.assertEqual(task.two()[0], 2.0)
    def test_loader_fail(self):
        self.taskfile.unlink()
        with self.assertRaises(IOError):
            loader.use(self.temptaskdir)
    def test_local_taskfile(self):
        local_taskfile = self.temptaskdir / 'taskfile.py'
        sample_taskfile.copy(local_taskfile)
        assert loader.which(self.temptaskdir).basename() == 'taskfile.py'
        task = loader.use(self.temptaskdir)
        self.assertEqual(task.two()[0], 2.0)
        local_taskfile.unlink()

