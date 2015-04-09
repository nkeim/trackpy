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
    """Implement a shadow "call stack" that tracks the chain of tasks.

    When an exception occurs, reports this information to the user.
    
    Instances of this class work as context managers.
    """
    def __init__(self, name, cwd=None):
        """Prepare a call stack entry for task 'name'
        
        'cwd' is the working directory. Only the cwd of the task
        that raised the exception is reported."""
        self.name = name
        self.cwd = cwd

    _tasker_stack = []  # Class attribute

    def __enter__(self):
        if self.name:
            self._tasker_stack.append(self.name)

    def __exit__(self, typ, val, tb):
        if tb is not None and len(self._tasker_stack):  # An exception occurred
            # Print some debug info above traceback
            sys.stderr.write('Tasker task(s): ' +
                    ' -> '.join(self._tasker_stack) + '\n')
            if self.cwd is not None:
                sys.stderr.write('Defined for "%s"\n' % self.cwd)
            if EDIT_TRACEBACKS:
                sys.stderr.write('NOTE: Some internal Tasker logic is not '
                        'shown in this traceback.\n'
                        'Set tasker.debug.EDIT_TRACEBACKS = False to '
                        'see everything.\n'
                        )
            sys.stderr.write('\n')
            # Print only one such message
            while self._tasker_stack: self._tasker_stack.pop()
        else:
            # Unwind our "call stack" naturally
            if self._tasker_stack: self._tasker_stack.pop()

