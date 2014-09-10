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

import sys

# Whether to hide implementation details in tracebacks
EDIT_TRACEBACKS = True


class tasker_traceback(object):
    def __init__(self, name, cwd=None):
        self.name = name
        self.cwd = cwd

    tasker_stack = []  # Class attribute

    def __enter__(self):
        if self.name:
            self.tasker_stack.append(self.name)

    def __exit__(self, typ, val, tb):
        if tb is not None and len(self.tasker_stack):  # An exception occurred
            # Print some debug info above traceback
            sys.stderr.write('Tasker task(s): ' +
                    ' -> '.join(self.tasker_stack) + '\n')
            if self.cwd is not None:
                sys.stderr.write('In %s\n' % self.cwd)
            # Print only one such message
            while self.tasker_stack: self.tasker_stack.pop()

