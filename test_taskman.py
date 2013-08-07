import os
import tempfile, unittest
from path import path
import pandas
from . import taskman, storage

basedir = path.getcwd()

def test_uniq():
    testdat = [0, 1, 1, 3, 7, 7, 7, 10, 10]
    from_uniq = taskman._uniq(testdat)
    from_sort = sorted(list(set(testdat)))
    if cmp(from_uniq, from_sort) != 0:
        raise AssertionError('%r != %r' % (from_uniq, from_sort))

class TaskerSubclass(taskman.Tasker):
    # Users will be encouraged to subclass Tasker, so we might as well do it here.
    # Here they could also define reader methods to read and process files in the
    # directory, outside of the pipeline (but perhaps invoking it).
    pass

class TestTaskman(unittest.TestCase):
    def enter_dir(self, dirname):
        """Set up a sample pipeline in a directory. 
        Tests depend on sample data defined here."""
        task = TaskerSubclass(dirname)
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
        def two(input):
            self.assertEqual(input, task.conf['one'])
            return task.conf['two'], {'twofloat': task.conf['two'], 
                    'onestr': input, 'name': task.conf['name']}
        @task([one, two], storage.Pandas('three.h5'))
        def three(ins):
            assert ins[0] == task.conf['one']
            twofloat = ins[1][1]['twofloat']
            assert twofloat == task.conf['two']
            return pandas.Series([twofloat,])
        @task({'three': three, 'td': 'three_dummy'}, 'four')
        def four(ins):
            assert ins['three'][0] == task.conf['two'] # First row of Series
            self.assertEqual(ins['td'].basename(), 'three_dummy')
            assert len(ins['td'].split()[0])
            assert not ins['td'].exists()
            (task.p / 'four').touch()
            return 'dummy'
        return task
    def setUp(self):
        os.chdir(basedir)
        self.testdir = path(tempfile.mkdtemp())
        #self.testdir = basedir / 'tmptest'
        if self.testdir.isdir():
            self.testdir.rmtree()
        self.testdir.makedirs_p()
        self.task = self.enter_dir(self.testdir)
    def tearDown(self):
        os.chdir(basedir)
        self.testdir.rmtree()
    def test_prepare_data(self):
        # Make sure we know what to do with TaskUnit and FileBase instances
        self.task.three()
        assert isinstance(self.task.one.outs[0], storage.FileBase)
        with self.task.p:
            self.assertEqual(self.task.one._prepare_data(self.task.one.outs_as_given), 
                    'one_str')
            self.assertEqual(taskman._nestmap(self.task.two._prepare_data, 
                self.task.one), 'one_str')
            # Returned data structure mirrors that of input
            self.assertIsInstance(taskman._nestmap(self.task.three._prepare_data, 
                self.task.three.ins_as_given), list)
    def test_2tasks(self):
        """Make sure up-to-date tasks do not rerun"""
        self.assertEqual(self.task.one(), 'one_str')
        self.assertEqual(self.task.two()[0], 2.0)
        self.assertEqual(self.task.one_count, 1) # two() should not re-run one()
    def test_pipeline(self):
        """Make sure full depth of task tree is traversed."""
        self.task.four() # Should be making one(), two(), etc.
        assert (self.testdir / 'one.json').exists()
    def test_load(self):
        self.task.two()
        self.task.two.load()
    def test_hash(self):
        assert len(set([self.task.one, self.task.one])) == 1
    def test_docstr(self):
        self.assertEqual(self.task.one.__doc__, 'Docstring: One')
    def test_clear(self):
        self.assertEqual(self.task.one(), 'one_str')
        self.task.two()
        self.assertEqual(self.task.one_count, 1) # two() should not re-run one()
        self.task.one.clear()
        self.assertEqual(self.task.one_count, 1)
        self.assertEqual(self.task.one(), 'one_str') # Re-runs one()
        self.assertEqual(self.task.one_count, 2)
    def test_report(self):
        self.assertSetEqual(set(self.task.four.report()), set(self.task.tasks.values()))
    def test_clearall(self):
        self.task.four()
        self.assertSetEqual(set(self.task.four.report()), set())
        self.task.clear()
        self.assertSetEqual(set(self.task.four.report()), set(self.task.tasks.values()))
    def test_twodirs(self):
        """Make sure separate task instances do not mix paths"""
        td2 = path(tempfile.mkdtemp())
        try:
            task2 = self.enter_dir(td2)
            self.assertNotEqual(self.task.two()[1]['name'], task2.two()[1]['name'])
            self.assertNotEqual(self.task.two()[1]['name'], task2.two()[1]['name'])
        finally:
            os.chdir(basedir)
            td2.rmtree()


