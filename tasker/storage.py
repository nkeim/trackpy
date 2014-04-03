import os
import json
import pandas
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
        try:
            hdf = pandas.HDFStore(self.filepath, 'r')
            r = hdf[self.key]
        finally:
            hdf.close()
        return r
    def save(self, data):
        self._mkdir()
        try:
            hdf = pandas.HDFStore(self.filepath, 'w')
            hdf[self.key] = data
        finally:
            hdf.close()

class JSON(FileBase):
    def read(self):
        return json.load(open(self.filepath, 'r'))
    def save(self, data):
        self._mkdir()
        json.dump(data, open(self.filepath, 'w'), indent=4, separators=(',', ': '))
