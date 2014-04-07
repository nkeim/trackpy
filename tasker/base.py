import os
from path import path

def _listify(arg):
    if isinstance(arg, (list, tuple)):
        return arg
    else:
        return [arg,]
def _toFiles(names):
    return [path(n) for n in names]

def isUpToDate(targets, deps=[], tartimes=[], deptimes=[], require=False):
    """Checks the modification times of files to see if targets need remaking.

    Each argument can be a list, or a single item.
    'targets' and 'deps' are lists of files. 
    'tartimes' and 'deptimes' are lists of additional mtimes to consider.

    Checks that all targets exist (even if no dependencies are given).
    If 'require' is set, raises an exception on missing dependencies.
    """
    targetFiles = _toFiles(_listify(targets))
    depFiles = _toFiles(_listify(deps))
    if not all([tf.exists() for tf in targetFiles]): return False
    if require and not all([df.exists() for df in depFiles]):
        raise IOError('Missing dependencies: %s' % \
                str([df for df in depFiles if not df.exists()]))
    targetFileTimes = [tar.mtime for tar in targetFiles]
    depFileTimes = [dep.mtime for dep in depFiles if dep.exists()]
    auxTargetTimes = _listify(tartimes)
    auxDepTimes = _listify(deptimes)
    if len(depFileTimes + auxDepTimes):
        return max(depFileTimes + auxDepTimes) <= min(targetFileTimes + auxTargetTimes)
    else:
        # No dependencies
        return True # Missing target files were already handled above.
class cachedprop(property):
    'Convert a method into a cached attribute'
    def __init__(self, method, doc=None):
        private = '_' + method.__name__
        def fget(s):
            try:
                return getattr(s, private)
            except AttributeError:
                value = method(s)
                setattr(s, private, value)
                return value
        def fdel(s):
            if private in s.__dict__:
                del s.__dict__[private]
        super(cachedprop, self).__init__(fget, fdel=fdel, doc=doc)
    @staticmethod
    def reset(self):
        cls = self.__class__
        for name in dir(cls):
            attr = getattr(cls, name)
            if isinstance(attr, cachedprop):
                delattr(self, name)
class DirBase(object):
    """Basis for directory-based data accessors.
    
    Attributes:
        'path' is what's passed to the constructor.
        'p' is what should be used to find enclosed files (instance of path.path).
        'opts' may be set through the constructor.

    Can also be used as a chdir() context manager.
    """
    def __init__(self, loc='.'):
        self.p = path(loc).abspath()
        self._chdir_stack = []
    name = property(lambda self: self.p.basename(),
            doc='Name of this directory')
    parentname = property(lambda self: (self.p / '..').basename(),
            doc='Name of parent directory')
    reset = cachedprop.reset
    def __repr__(self):
        return '%s: %s' % (self.__class__.__name__, self.p)
    def __str__(self):
        return str(self.p)
    def __enter__(self):
        self._chdir_stack.append(os.getcwd())
        os.chdir(self.p)
    def __exit__(self, exc_type, exc_val, exc_tb):
        os.chdir(self._chdir_stack.pop(-1))

