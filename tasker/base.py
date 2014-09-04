#   Copyright 2014 Nathan C. Keim
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

import os
from path import path

class AttrDict(dict):
    """Dictionary that also exposes items as attributes."""
    def __init__(self, *args, **kw):
        super(AttrDict, self).__init__(*args, **kw)
        self.__dict__ = self


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

    name = property(lambda self: str(self.p.basename()),
            doc='Name of this directory')

    parentname = property(lambda self: str((self.p / '..').basename()),
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

