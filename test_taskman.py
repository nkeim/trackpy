import sys
import tempfile, shutil, unittest
from path import path
import pandas
from . import taskman, storage

def test_uniq():
    testdat = [0, 1, 1, 3, 7, 7, 7, 10, 10]
    from_uniq = taskman._uniq(testdat)
    from_sort = sorted(list(set(testdat)))
    if cmp(from_uniq, from_sort) != 0:
        raise AssertionError('%r != %r' % (from_uniq, from_sort))

class TaskManSubclass(taskman.TaskMan):
    # Users will be encouraged to subclass TaskMan, so we might as well do it here.
    # Here they could also define reader methods to read and process files in the
    # directory, outside of the pipeline (but perhaps invoking it).
    pass

class test_Taskman(unittest.TestCase):
    def enter_dir(self, dirname):
        # Set up a sample pipeline in a directory -- part of setup code
        task = TaskManSubclass(dirname)
        # Totally arbitrary configuration "scheme"
        task.conf = dict(one='one_str', two=2.0, name=task.name) 
        task.one_count = 0
        @task([], storage.JSON('one.json'))
        def one(ins):
            """Docstring: One"""
            task.one_count += 1
            assert len(ins) == 0
            return task.conf['one'] # Goes to JSON
        @task(one, [storage.JSON('two.json'), storage.JSON('2b.json')])
        def two(ins):
            """Docstring: Two"""
            self.assertEqual(ins[0], task.conf['one']) # Index by position
            assert ins['one'] == task.conf['one'] # Index by filename
            assert ins['one.json'] == task.conf['one'] # Index by filename
            return task.conf['two'], {'twofloat': task.conf['two'], 'onestr': ins['one']}
        @task([one, two], storage.Pandas('three.h5'))
        def three(ins):
            assert task.one_count == 1 # two() should not re-run one()
            assert ins[0]['one'] == task.conf['one']
            twofloat = ins[1]['2b']['twofloat']
            assert twofloat == task.conf['two']
            return pandas.Series([twofloat,])
        @task((three, 'three_dummy'), 'four')
        def four(ins):
            assert ins[0]['three'][0] == task.conf['two'] # First row of Series
            self.assertEqual(ins[1].basename(), 'three_dummy')
            assert len(ins[1].split()[0])
            assert not ins[1].exists()
            return 'dummy'
        return task
    def setUp(self):
        self.testdir = path(tempfile.mkdtemp())
        self.task = self.enter_dir(self.testdir)
    def test_prepare_data(self):
        # Make sure we know what to do with TaskUnit and FileBase instances
        self.task.one()
        assert isinstance(self.task.one.outs[0], storage.FileBase)
        with self.task.p:
            self.assertEqual(self.task.one._prepare_data(self.task.one.outs[0]), 'one_str')
            self.assertIsInstance(self.task.two._prepare_data(self.task.one), dict)
            self.assertIsInstance(self.task.two._prepare_data(self.task.two.ins_as_given),
                    dict)
            # Returned data structure mirrors that of input
            self.assertIsInstance(self.task.two._prepare_data(self.task.two.ins), list)
    def test_2tasks(self):
        self.assertEqual(self.task.one(), 'one_str')
        self.assertEqual(self.task.two()[0], 2.0)
    def test_pipeline(self):
        self.task.four() # Should be making one(), two(), etc.
        assert (self.testdir / 'one.json').exists()
    @unittest.skip
    def test_load(self):
        self.task.two()
        self.task.two.load()
    def test_hash(self):
        assert len(set([self.task.one, self.task.one])) == 1
    def test_docstr(self):
        self.assertEqual(self.task.one.__doc__, 'Docstring: One')
    def tearDown(self):
        shutil.rmtree(self.testdir)


