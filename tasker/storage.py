import os
import json
import pandas
from path import path

class FileBase(object):
    def __init__(self, filename, parentdir='.'):
        """Specify location of file. If 'filename' is absolute, 'parentdir' is ignored."""
        self.filename = path(filename)
        self.set_parentdir(parentdir)
    def set_parentdir(self, parentdir):
        self.parentdir = path(parentdir)
        if self.filename.isabs():
            self.filepath = self.filename
        else:
            self.filepath = (self.parentdir / self.filename).normpath().abspath()
    def read(self): pass # Uses self.filepath
    def save(self, data): pass # Uses self.filepath
    def __call__(self):
        return self.read()

class Pandas(FileBase):
    def read(self):
        try:
            hdf = pandas.HDFStore(self.filepath, 'r')
            r = hdf[os.path.splitext(os.path.basename(self.filepath))[0]]
        finally:
            hdf.close()
        return r
    def save(self, data):
        try:
            hdf = pandas.HDFStore(self.filepath, 'w')
            hdf[os.path.splitext(os.path.basename(self.filepath))[0]] = data
        finally:
            hdf.close()

class JSON(FileBase):
    def read(self):
        return json.load(open(self.filepath, 'r'))
    def save(self, data):
        json.dump(data, open(self.filepath, 'w'), indent=4, separators=(',', ': '))
