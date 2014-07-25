import os
import tempfile, unittest
from path import path
import pandas
from tasker import task, storage, progress

basedir = path.getcwd()

def test_uniq():
    testdat = [0, 1, 1, 3, 7, 7, 7, 10, 10]
    from_uniq = task._uniq(testdat)
    from_sort = sorted(list(set(testdat)))
    if cmp(from_uniq, from_sort) != 0:
        raise AssertionError('%r != %r' % (from_uniq, from_sort))

class TaskerSubclass(task.Tasker):
    # Users will be encouraged to subclass Tasker, so we might as well do it here.
    # Here they could also define reader methods to read and process files in the
    # directory, outside of the pipeline (but perhaps invoking it).
    pass

class TestTask(unittest.TestCase):
    def enter_dir(self, dirname):
        """Set up a sample pipeline in a directory. 
        Tests depend on sample data defined here."""
        task = TaskerSubclass(dirname)
        # Totally arbitrary configuration "scheme"
        task.conf = dict(one='one_str', two=2.0, name=task.name) 
        task.one_count = 0
        @task.create_task([], storage.JSON('one.json'))
        def one(tsk, ins):
            """Docstring: One"""
            task.one_count += 1
            assert len(ins) == 0
            return task.conf['one'] # Goes to JSON
        @task.create_task(one, [storage.JSON('two.json'), storage.JSON('2b.json')])
        def two(tsk, input):
            self.assertEqual(input, task.conf['one'])
            return task.conf['two'], {'twofloat': task.conf['two'], 
                    'onestr': input, 'name': task.conf['name']}
        @task.create_task([one, two], storage.Pandas('three.h5'))
        def three(tsk, ins):
            assert ins[0] == task.conf['one']
            twofloat = ins[1][1]['twofloat']
            assert twofloat == task.conf['two']
            return pandas.Series([twofloat,])
        @task.create_task({'three': three, 'td': 'three_dummy'}, 'four')
        def four(tsk, ins):
            assert ins['three'][0] == task.conf['two'] # First row of Series
            self.assertEqual(ins['td'].basename(), 'three_dummy')
            assert len(ins['td'].split()[0])
            (task.p / 'four').touch()
            return 'dummy'
        return task
    def setUp(self):
        os.chdir(basedir)
        self.testdir = path(tempfile.mkdtemp())
        self.task = self.enter_dir(self.testdir)
    def tearDown(self):
        os.chdir(basedir)
        self.testdir.rmtree()
    def test_name(self):
        assert self.task.one.__name__ == 'one'
    def test_prepare_data(self):
        # Make sure we know what to do with TaskUnit and FileBase instances
        self.task.three()
        assert isinstance(self.task.one.outs[0], storage.FileBase)
        with self.task:
            self.assertEqual(self.task.one._prepare_data(self.task.one._outs_as_given), 
                    'one_str')
            self.assertEqual(task._nestmap(self.task.two._prepare_data,
                self.task.one), 'one_str')
            # Returned data structure mirrors that of input
            self.assertIsInstance(task._nestmap(self.task.three._prepare_data,
                self.task.three._ins_as_given), list)
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
        assert self.task.one() == 'one_str'
        self.task.two()
        assert self.task.one_count == 1  # two() should not re-run one()
        self.task.one.clear()
        assert self.task.one_count == 1
        assert self.task.one() == 'one_str'  # Re-runs one()
        assert self.task.one_count == 2
    def test_report(self):
        self.assertSetEqual(set(self.task.four.report()),
                            set([self.task.one, self.task.two,
                                 self.task.three, self.task.four]))
    def test_clearall(self):
        self.task.four()
        self.assertSetEqual(set(self.task.four.report()), set())
        self.task.clear()
        self.assertSetEqual(set(self.task.four.report()),
                            set([self.task.one, self.task.two,
                                 self.task.three, self.task.four]))
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

    def test_context(self):
        """Text chdir context behavior"""
        assert path('.').abspath().basename() != self.testdir.abspath().basename()
        with self.task:
            assert path('.').abspath().basename() == self.testdir.abspath().basename()

    def test_progress(self):
        @self.task.create_task([], ['progress.tag'])
        def exercise_progress(tsk, ins):
            mon = progress.Monitor([self.task.p])
            def readprog():
                assert (self.task.p / progress.DEFAULT_STATUS_FILE).exists()
                stat = mon.get_statuses()
                return stat.ix[0].to_dict()
            assert readprog()['status'] == 'working'
            assert 'elapsed_time' in readprog()
            for i in tsk.progress.tally(range(10)):
                assert type(i) is int
                stat = readprog()
                assert stat['status'] == 'working'
                assert stat['current'] == i + 1
                assert 'time_per' in stat
            for i in tsk.progress.tally(range(10), 10):
                assert type(i) is int
                stat = readprog()
                assert stat['total'] == 10
                assert stat['current'] == i + 1
                assert 'time_per' in stat
            for i in range(10):
                tsk.progress.working(i)
                stat = readprog()
                assert stat['status'] == 'working'
                assert stat['current'] == i + 1
                assert 'time_per' not in stat
            for i in range(10):
                tsk.progress.working(i, 10)
                stat = readprog()
                assert stat['current'] == i + 1
            for i in range(10):
                tsk.progress.working(i, 10, {'hi': 1.5})
                stat = readprog()
                assert stat['current'] == i + 1
                assert stat['hi'] == 1.5
        self.task.exercise_progress()

    def test_locking(self):
        didrun = []
        @self.task.create_task([], ['dummy.tag'])
        def try_locking(tsk, ins):
            didrun.append(True)
            assert self.task.is_working()
            assert self.task.is_working(task='try_locking')
            assert not self.task.is_working(task='something_else')
            self.task.one.run() # Should be OK
            self.assertRaises(task.LockException, try_locking)
        assert not self.task.is_working()
        self.task.try_locking()
        assert didrun
        assert not self.task.is_working()
        self.task.unlock() # Just to try it
        self.task.unlock() # Just to try it again
        assert not self.task.is_working()

    def test_menu(self):
        self.task.menu()
        self.assertEqual(self.task.one(), 'one_str')
        # Check display of completed task
        self.task.menu()


class TestNewStyleTasks(TestTask):
    def enter_dir(self, dirname):
        """Set up a sample pipeline in a directory.
        Tests depend on sample data defined here."""
        task = TaskerSubclass(dirname)
        # Totally arbitrary configuration "scheme"
        task.conf = dict(one='one_str', two=2.0, name=task.name)
        task.one_count = 0
        @task.stores(storage.JSON('one.json'))
        def one(tsk):
            """Docstring: One"""
            task.one_count += 1
            return task.conf['one'] # Goes to JSON
        @task.stores(storage.JSON('two.json'), storage.JSON('2b.json'))
        def two(tsk, one=one):
            self.assertEqual(one, task.conf['one'])
            return [task.conf['two'],
                        {'twofloat': task.conf['two'],
                            'onestr': one, 'name': task.conf['name']}]
        @task.stores(storage.Pandas('three.h5')) # Positional form
        def three(tsk, one=one, two=two):
            assert one == task.conf['one']
            twofloat = two[1]['twofloat']
            assert twofloat == task.conf['two']
            return pandas.Series([twofloat,])
        @task.stores('four')
        def four(tsk, three=three, td='three_dummy'):
            assert three[0] == task.conf['two'] # First row of Series
            self.assertEqual(td.basename(), 'three_dummy')
            assert len(td.split()[0])
            (task.p / 'four').touch()
            return 'dummy'

        # No storage
        @task
        def doesnt_store(tsk, three=task.three):
            return three[0]
        @task.computes # Alternate syntax
        def doesnt_store2(tsk, three_val=task.doesnt_store):
            return three_val

        # Storage with no-store gap
        @task.stores(storage.JSON('gapped.json'))
        def gapped(tsk, three_val=task.doesnt_store):
            return three_val

        return task

    def test_nostore_task(self):
        assert not (self.task.p / 'three.h5').exists()
        assert self.task.doesnt_store() == self.task.conf['two']
        assert (self.task.p / 'three.h5').exists()
        assert self.task.doesnt_store2() == self.task.conf['two']

    def test_nostore_gap(self):
        """Check whether updates to upstream tasks propagate through
        a no-store task.
        """
        assert self.task.gapped() == self.task.conf['two']
        self.task.conf['two'] = 100
        self.task.two.clear()
        assert self.task.gapped() == self.task.conf['two']

    def test_nostore_current(self):
        """Check whether a no-store task correctly reports its status.
        
        (Which should reflect the status of upstream tasks.)"""
        assert not self.task.three.is_current()
        assert not self.task.doesnt_store.is_current()
        # Run task
        assert self.task.doesnt_store() == self.task.conf['two']
        assert self.task.three.is_current()
        assert self.task.doesnt_store.is_current()
        assert not self.task.gapped.is_current()
        assert len(self.task.gapped.report()) == 1
        assert self.task.gapped() == self.task.conf['two']
        assert self.task.gapped.is_current()
        assert self.task.doesnt_store.is_current()

    def test_nostore_clear(self):
        """clear() on a no-storing task should raise an error."""
        self.task.gapped()
        assert self.task.gapped.is_current()
        assert self.task.doesnt_store.is_current()
        self.assertRaises(NotImplementedError,
                          self.task.doesnt_store.clear)

    def test_prepare_data(self):
        """Make sure we know what to do with TaskUnit and FileBase instances"""
        self.task.three()
        assert isinstance(self.task.one.outs[0], storage.FileBase)
        with self.task:
            self.assertEqual(self.task.one._prepare_data(self.task.one._outs_as_given),
                             'one_str')
            self.assertEqual(task._nestmap(self.task.two._prepare_data,
                                           self.task.one), 'one_str')
            # Returned data structure mirrors that of input
            self.assertIsInstance(task._nestmap(self.task.three._prepare_data,
                                                self.task.three._ins_as_given), dict)

    def test_file_dependency(self):
        """A FileBase object is specified instead of a task."""
        prebound = storage.JSON(self.testdir / 'prebound.json')
        @self.task
        def use_files(tsk, pb=prebound,
                         runtime_bound=storage.JSON('runtime_bound.json')):
            return True
        input_filenames = [fn.basename() for fn in use_files.input_files]
        assert len(input_filenames) == 2
        assert 'prebound.json' in input_filenames
        assert 'runtime_bound.json' in input_filenames

    def test_output_sequence(self):
        """Check for bug in which a task returning a sequence couldn't
        store to a single file."""
        @self.task.stores('a_list.json')
        def a_list(tsk):
            return [1, 2, 3]
        a_list()

    def test_traceback(self):
        @self.task
        def raises(tsk):
            raise RuntimeError()
        @self.task
        def deps_on_raises(tsk, r=raises):
            pass
        try:
            deps_on_raises()
        except RuntimeError:
            pass

