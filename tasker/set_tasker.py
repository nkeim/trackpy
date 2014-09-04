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

# FIXME: No test coverage

from path import path

from .base import cachedprop
from .task import Tasker
from .storage import JSON
from .loader import use

class SetTasker(Tasker):
    """Tasker with convenient methods for handling taskers in subdirectories"""
    def __init__(self, *args, **kw):
        super(SetTasker, self).__init__(*args, **kw)

    def use(self, dirname, **kw):
        """use() a directory with path relative to this one's"""
        return use(self.p / dirname, **kw)

    @cachedprop
    def groups(self):
        """Dictionary of lists of relative paths to subdirectories."""
        return JSON(self.p / 'groups.json').read()

    def group_paths(self, groupname):
        """Get absolute paths of a named group of subdirectories"""
        return [path(t).abspath() for t in self.groups[groupname]]

    def grp(self, groupname, **kw):
        """use() each directory in a named group"""
        return [use(self.p / dirname, **kw) for dirname in self.groups[groupname]]

    @cachedprop
    def aliases(self):
        """Dictionary of aliases to relative paths of subdirectories."""
        return JSON(self.p / 'aliases.json').read()

    def alias_path(self, aliasname):
        """Get absolute path of an aliased subdirectory"""
        return [path(t).abspath() for t in self.aliases[aliasname]]

    def al(self, aliasname, **kw):
        """use() an aliased subdirectory"""
        return use(self.p / self.aliases[aliasname], **kw)

