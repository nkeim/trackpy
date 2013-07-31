import os
import json
import pandas
from path import path

class FileBase(object):
    def __init__(self, filename):
        self.filename = path(filename)
    def read(self): pass
    def save(self, data): pass
class Pandas(FileBase):
    def read(self):
        try:
            hdf = pandas.HDFStore(self.filename, 'r')
            r = hdf[os.path.splitext(os.path.basename(self.filename))[0]]
        finally:
            hdf.close()
        return r
    def save(self, data):
        try:
            hdf = pandas.HDFStore(self.filename, 'w')
            hdf[os.path.splitext(os.path.basename(self.filename))[0]] = data
        finally:
            hdf.close()
class JSON(FileBase):
    def read(self):
        return json.load(open(self.filename, 'r'))
    def save(self, data):
        json.dump(data, open(self.filename, 'w'), indent=4, separators=(',', ': '))
