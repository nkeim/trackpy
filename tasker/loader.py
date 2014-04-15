import imp, hashlib, sys
from path import path

def use(directory, taskfile='taskfile.py', taskfile_sub='taskfile_sub.py', **kw):
    """Return an object for 'directory' that gives access to its tasks.

    This is an instance of Tasker (or a subclass thereof) that is defined and
    set up in the use() function, contained in 'taskfile' (default "taskfile.py").

    If 'taskfile' is not in 'directory', searches upward from 'directory' to find
    a file named 'taskfile_sub' (default "taskfile_sub.py"). 'taskfile' and
    'taskfile_sub' can be glob patterns. They can also be absolute paths,
    in which case no search is performed.
    
    The use() that is found is called with 'directory' and any additional keyword 
    arguments.
    """
    tfmod = taskmod(directory, taskfile=taskfile)
    try:
        return tfmod.use(directory, **kw)
    except:
        sys.stderr.write('Problem calling use() for %s\n' % directory)
        raise

def taskmod(directory, taskfile='taskfile.py', taskfile_sub='taskfile_sub.py'):
    """Import the module that defines tasks (and anything else) for 'directory'.

    If 'taskfile' is not in 'directory', searches upward from 'directory' to find
    a file named 'taskfile_sub' (default "taskfile_sub.py"). 'taskfile' and
    'taskfile_sub' can be glob patterns. They can also be absolute paths,
    in which case no search is performed.
    """
    tf = which(directory, taskfile=taskfile, taskfile_sub=taskfile_sub)
    dirpath = path(directory).abspath()
    # Assign a unique module name, to avoid trouble
    tf_hash = hashlib.sha1(str(tf) + str(dirpath)).hexdigest()[:10]
    return imp.load_source('_taskfile_' + tf_hash, tf)

def which(directory, taskfile='taskfile.py', taskfile_sub='taskfile_sub.py'):
    """Returns path to the Python file defining tasker tasks for 'directory'.

    'taskfile' is the glob pattern to match with filenames, or an absolute path.
    """
    if path(taskfile).isabs(): 
        return path(taskfile)
    dirpath = path(directory).abspath()
    tfs_local = dirpath.glob(taskfile)
    if len(tfs_local) == 0:
        tf = _find_upward(dirpath.parent, taskfile_sub)
        if tf is None:
            raise IOError('Task definitions not found in %s or any directory above it' \
                    % (dirpath))
        else:
            return tf
    elif len(tfs_local) == 1:
        return tfs_local[0]
    else:
        raise IOError('Found multiple matches for %s in %s: %r' \
                      % (taskfile, dirpath, tfs_local))

def _find_upward(dirpath, pattern):
    """Search for glob pattern 'pattern' going up the directory tree.
    
    Returns path of first file found, or None if the top of the tree was reached.
    If 'pattern' is an absolute pathname, does not perform a search.
    
    Raises exception if multiple matching files are found in a single directory."""
    if path(pattern).isabs():
        return path(pattern)
    p = dirpath
    while p != p.parent: # Stop when we can't go up any further (at root path)
        tfs = p.glob(pattern)
        if len(tfs) == 0:
            p = p.parent
        elif len(tfs) == 1:
            return tfs[0]
        else:
            raise IOError('Found multiple matches for %s in %s: %r' % (pattern, p, tfs))
    else:
        return None

