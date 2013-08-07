import os
from path import path

from .base import isUpToDate, DirBase
from .storage import FileBase

def _listify(arg):
    if isinstance(arg, (list, tuple)):
        return arg
    else:
        return [arg,]
def _toFiles(names):
    return [path(n) for n in names]
def _uniq(l):
    """Removes duplicates from a list, preserving order."""
    r = []
    for item in l:
        if item not in r: r.append(item)
    return r
class TaskUnit(object):
    def __init__(self, func, ins, outs, taskman):
        """Create a task.
        'func' is a function which turns 'ins' into 'outs'.
        'ins' is some data structure (preferably list or dict) populated
            with filenames, FileBase instances, or other TaskUnit instances.
        'outs' is a filename or FileBase instance, or a list thereof.
        'taskman' is a parent TaskMan instance, which presently serves to 
            set the working directory for this task.

        When 'func' is called, it takes a single argument which is structured
        like 'ins', but with FileBase instances replaced with those files' contents,
        and with each TaskUnit instance replaced with a dict that contains that TaskUnit's
        outputs, accessed through enumeration, basename, or full filename.
        'func' returns a list which corresponds to the elements of 'outs'. Elements
        not corresponding to a FileBase instance are ignored.

        In general, any filename can be either relative to the task's working directory.
        or absolute. 'func' will always be executed in the task's working directory.

        When a user asks for the outputs of 'func', they are loaded from disk and 
        returned as a dict, which is indexed just like in the argument to 'func'.
        """
        self.func = func
        self.__name__ = func.__name__
        self.__doc__ = func.__doc__
        self.taskman = taskman
        self.p = self.taskman.p # This directory will always be my working dir
        self.outs_as_given = outs
        self.outs = _listify(outs)
        self.output_files = map(self._get_filename, self.outs)
        self.ins_as_given = ins
        self.ins = _listify(ins)
        self.input_files, self.input_tasks = self._flatten_dependencies()
    # Deal with arbitrary user specification of task inputs
    def _get_filename(self, fileobj):
        """Returns absolute path to a file, from a string or FileBase."""
        if isinstance(fileobj, FileBase): 
            fn = path(fileobj.filename)
        else:
            fn = path(fileobj)
        if fn.isabs(): return fn
        else: return (self.p / fn).abspath()
    def _prepare_data(self, data_part):
        """Returns input data for self.func(), mirroring the structure of 'ins'.
        
        This needs to be run in my working dir."""
        if isinstance(data_part, (str, path)): # Literal filename
            return path(data_part).abspath()
        elif isinstance(data_part, FileBase): # Formatted file
            return data_part.read()
        elif isinstance(data_part, TaskUnit): # Another task
            d = {} # To be populated with various aliases for the data
            for i, indatum in enumerate(data_part.outs):
                infile = data_part._get_filename(indatum)
                data = self._prepare_data(indatum) # goes to _load_file()
                d.update({k: data for k in (os.path.basename(infile), \
                        os.path.splitext(os.path.basename(infile))[0], i)})
            return d
        # Cases to allow nesting
        elif isinstance(data_part, dict):
            return {k: self._prepare_data(data_part[k]) for k in data_part}
        elif isinstance(data_part, (list, tuple)):
            return [self._prepare_data(v) for v in data_part]
        else:
            raise ValueError('%r does not specify a valid input source' % data_part)
    def _flatten_dependencies(self):
        """Returns a list of files and a list of tasks.
        Files are returned as absolute paths."""
        deps = self._flatten_dependencies_recurse(self.ins)
        files, tasks = set(), set()
        for v in deps:
            if isinstance(v, TaskUnit): tasks.add(v)
        # File paths should be relative to my working dir.
        return [f if f.isabs() else (self.p / f).abspath() \
                for f in _toFiles(list(files))], \
                list(tasks)
    def _flatten_dependencies_recurse(self, in_part):
        """Walks the 'self.ins' data structure in a manner similar to _prepare_inputs()"""
        flatten = lambda seq: reduce((lambda a,b: a+b), seq, [])
        if isinstance(in_part, (str, FileBase, TaskUnit)):
            return [in_part,]
        elif isinstance(in_part, dict):
            return flatten([self._flatten_dependencies_recurse(v) \
                    for v in in_part.values()])
        elif isinstance(in_part, (list, tuple)):
            return flatten([self._flatten_dependencies_recurse(v) for v in in_part])
        else:
            raise ValueError('Part of input specification could not be handled: %s' \
                    % repr(in_part))
    def __repr__(self):
        return 'TaskUnit: ' + self.func.__name__
    def __hash__(self):
        return id(self) # So we can put these in sets
    # Public interface
    def __call__(self):
        """Update outputs if necessary and read from disk. 
        
        See load() for details of return value."""
        self.sync()
        return self.load()
    def load(self):
        """Return representation of task output on disk.

        Structure is the same as in the task definition --- a single
        value or a sequence of values.

        Where a FileBase instance was used (i.e. a recognized file format like JSON()),
        the contents of the file are returned. Otherwise, the path to the file is
        returned.
        """
        # We return a dict in which the values are referred to by various names,
        # the wame way we pass input data to self.func() itself.
        with self.p: return self._prepare_data(self.outs_as_given)
    def run(self):
        """Execute task (always) and write output files in recognized formats."""
        with self.p:
            outdata = _listify(self.func(self._prepare_data(self.ins_as_given)))
            assert len(outdata) == len(self.outs)
            for of, od in zip(self.outs, outdata):
                if isinstance(of, FileBase): of.save(od)
    def sync(self):
        """run() task if required to keep outputs up to date."""
        for it in self.input_tasks: it.sync() # Really, really inefficient
        if not isUpToDate(self.output_files, self.input_files): self.run()
    def isUpToDate(self):
        """True if this task and all its dependencies are current."""
        return isUpToDate(self.output_files, self.input_files) and \
                all([it.isUpToDate() for it in self.input_tasks])
    def report(self):
        """List tasks that would have to be run to update outputs."""
        return _uniq(self._report_recurse())
    def _report_recurse(self):
        """Returns list of sub-tasks that are out of date."""
        r = reduce(lambda a,b: a+b, [[],] + [t._report_recurse() for t in self.input_tasks])
        if not isUpToDate(self.output_files, self.input_files): r.append(self)
        return r
    def clear(self):
        """Delete this task's output files and directories."""
        for f in self.output_files():
            if f.isdir(): f.rmtree()
            elif f.isfile(): f.unlink()
class TaskMan(DirBase):
    """Object to set up tasks within a single directory.
    
    Initialize with the directory name."""
    def __init__(self, dirname):
        super(TaskMan, self).__init__(dirname)
        self.tasks = {}
        self.conf = {}
    def __call__(self, ins=None, outs=None):
        """Returns a decorator that turns the function into a task 
        with registered inputs and outputs."""
        if ins is None: ins = []
        if outs is None: outs = []
        def mktask(func): 
            t = TaskUnit(func, ins, outs, self)
            self.tasks[t.__name__] = t
            setattr(self, t.__name__, t)
            return t
        return mktask
    def which(self, filename):
        """Return the TaskUnit instance that is responsible for 'filename',
        which can be relative to this instance's working dir, or absolute."""
        fnp = path.filename
        if fnp.isabs: fnp_abs = fnp
        else: fnp_abs = (self.p / fnp).abspath()
        for t in self.tasks:
            if fnp_abs in t.output_files:
                return t
    def clear(self):
        """Removes all output files of all tasks."""
        for t in self.tasks:
            t.clear()

