import os
import cPickle
import json
from path import path

class FileBase(object):
    def __init__(self, filename, parentdir='.'):
        """Specify location of file. If 'filename' is absolute, 'parentdir' is ignored."""
        self.filename = path(filename).normpath()
        self.set_parentdir(parentdir)
    def set_parentdir(self, parentdir):
        self.parentdir = path(parentdir)
        if self.filename.isabs():
            self.filepath = self.filename
        else:
            self.filepath = (self.parentdir / self.filename).normpath().abspath()
    def read(self): pass # Uses self.filepath
    def save(self, data): pass # Uses self.filepath
    def _mkdir(self):
        """Makes directory that self.filepath goes in, if it does
        not exist.
        """
        self.filepath.dirname().makedirs_p()
    def __call__(self):
        return self.read()

class Pandas(FileBase):
    """Store a Pandas data object."""
    def __init__(self, filename, parentdir='.', key=None):
        """Specify location of file. 
        
        If 'filename' is absolute, 'parentdir' is ignored.

        'key' gives the name of the Pandas object within the file.
            Defaults to the filename (minus extension).
        """
        super(Pandas, self).__init__(filename, parentdir=parentdir)
        if key is None:
            self.key = os.path.splitext(os.path.basename(self.filepath))[0]
        else:
            self.key = key

    def read(self):
        import pandas
        hdf = pandas.HDFStore(self.filepath, 'r')
        try:
            r = hdf[self.key]
        finally:
            hdf.close()
        return r

    def save(self, data):
        import pandas
        self._mkdir()
        hdf = pandas.HDFStore(self.filepath, 'w')
        try:
            hdf[self.key] = data
        finally:
            hdf.close()

class JSON(FileBase):
    """Store basic Python types (dicts, lists, floats, ints, etc.)
    in a human-readable text format.
    """
    def read(self):
        with open(self.filepath, 'r') as f:
            return json.load(f)

    def save(self, data):
        self._mkdir()
        with open(self.filepath, 'w') as f:
            json.dump(data, f, indent=4, separators=(',', ': '))

class Pickle(FileBase):
    """Store most any Python object in a Python-only binary format.
    """
    def read(self):
        with open(self.filepath, 'r') as f:
            return cPickle.load(f)

    def save(self, data):
        self._mkdir()
        with open(self.filepath, 'w') as f:
            cPickle.dump(data, f)
