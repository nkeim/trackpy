import imp, hashlib, sys
from path import path

def use(directory, taskfile='taskfile.py', **kw):
    """Return an object for 'directory' that gives access to its tasks.

    This is an instance of Tasker (or a subclass thereof) that is defined and
    set up in the use() function, contained in 'taskfile' (default "taskfile.py").

    Searches upward from 'directory' to find the relevant 'taskfile'. 'taskfile'
    can be a glob pattern or an absolute path, in which case no search is performed.
    
    The use() that is found is called with 'directory' and any additional keyword 
    arguments.
    """
    tfmod = taskmod(directory, taskfile=taskfile)
    try:
        return tfmod.use(directory, **kw)
    except:
        sys.stderr.write('Problem calling use() for %s\n' % directory)
        raise

def taskmod(directory, taskfile='taskfile.py'):
    """Import the module that defines tasks (and anything else) for 'directory'.

    Searches upward from 'directory' to find the relevant 'taskfile'. 'taskfile'
    can be a glob pattern or an absolute path, in which case no search is performed.
    """
    tf = which(directory, taskfile=taskfile)
    # Assign a unique module name, to avoid trouble
    tf_hash = hashlib.sha1(str(tf) + str(directory)).hexdigest()[:10]
    return imp.load_source('_taskfile_' + tf_hash, tf)

def which(directory, taskfile='taskfile.py'):
    """Returns path to the Python file defining tasker tasks for 'directory'.

    'taskfile' is the glob pattern to match with filenames, or an absolute path.
    """
    if path(taskfile).isabs(): 
        return path(taskfile)
    dirpath = path(directory).abspath()
    tf = _find_upward(dirpath, taskfile)
    if tf is None:
        raise IOError('Task definitions ("%s") not found in %s or any directory above it' \
                % (taskfile, dirpath))
    return tf

def _find_upward(dirpath, pattern):
    """Search for glob pattern 'pattern' going up the directory tree.
    
    Returns path of first file found, or None if the top of the tree was reached.
    
    Raises exception if multiple matching files are found in a single directory."""
    p = dirpath
    while p != p.parent: # Stop when we reach the root path
        tfs = p.glob(pattern)
        if len(tfs) == 0:
            p = p.parent
        elif len(tfs) == 1:
            return tfs[0]
        else:
            raise IOError('Found multiple matches for %s in %s: %r' % (pattern, p, tfs))
    else:
        return None

