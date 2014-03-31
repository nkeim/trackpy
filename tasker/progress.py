"""Support functions for communicating worker status."""

import os, json, time, datetime
import signal
import numpy as np
import pandas

class Monitor(object):
    """Monitor workers in many directories."""
    def __init__(self, dirs, filename='taskerstatus.json'):
        self.dirs = dirs
        self.filename = filename
    def get_statuses(self):
        """Returns DataFrame of status info for a list of filenames"""
        info = []
        for d in self.dirs:
            try:
                sfn = os.path.join(d, self.filename)
                sf = open(sfn, 'r')
                sfinfo = json.load(sf)
                sfinfo['dir'] = str(d)
                sfinfo['since_update'] = _format_td(datetime.timedelta(0,
                    time.time() - os.path.getmtime(sfn)))
            except IOError:
                sfinfo = {'dir': str(d), 'status': '?'}
            info.append(sfinfo)
        return pandas.DataFrame(info)
    def show(self, custom_columns=None):
        """Presents formatted status info.

        'custom_columns' : list of column names to display

        Returns a DataFrame, which should display nicely.
        """
        df = self.get_statuses()
        columns = ['dir', 'task', 'since_update', 'status',
                   'elapsed_time']
        if custom_columns is None:
            columns.extend([ 'current', 'total', 'time_left',])
        else:
            columns.extend(custom_columns)
        for cn in columns:
            if cn not in df:
                df[cn] = ''
        return df[columns]
    def watch(self, interval=5, custom_columns=None):
        """Displays regularly refreshed status_board() in IPython.

        'interval' : update interval, in seconds

        Catches a KeyboardInterrupt (e.g. Ctrl-m Ctrl-i in the notebook)
        and exits cleanly.
        """
        import IPython.display
        try:
            while True:
                sb = self.show(custom_columns=custom_columns)
                IPython.display.clear_output()
                IPython.display.display_html(sb.to_html(na_rep=''), raw=True)
                time.sleep(interval)
        except KeyboardInterrupt:
            IPython.display.clear_output()
            IPython.display.display_html(sb.to_html(na_rep=''), raw=True)
            print 'Last update: ' + datetime.datetime.now().strftime('%c')
            return
    def abort(self, indices=None):
        """Interrupts tasks running LOCALLY on this computer.

        indices : Optional sequence of indices to the DataFrame displayed by show()
        """
        if indices is not None:
            pids = self.get_statuses().ix[list(indices)].pid
        else:
            pids = self.get_statuses().pid
        for pid in pids:
            try:
                pid = int(pid)
            except ValueError:
                pass
            else:
                os.kill(pid, signal.SIGINT)

class StatusFile(object):
    """JSON-formatted file for status of a long-running computation"""
    def __init__(self, persistent_info=None, filename='taskerstatus.json'):
        """'persistent_info' is a dict that will be included
        in every update.
        """
        self.filename = filename
        if persistent_info is None:
            self.persistent_info = {}
        else:
            self.persistent_info = persistent_info.copy()
    def update(self, newinfo):
        """Write status file with 'newinfo', including persistent information."""
        tmpname = self.filename + '._tmp'
        tmpfile = open(tmpname, 'w')
        info = self.persistent_info.copy()
        info.update(newinfo)
        json.dump(info, tmpfile, indent=4, separators=(',', ': '))
        tmpfile.close()
        if os.name == 'nt':
            os.unlink(self.filename) # Windows doesn't allow overwriting existing file
        os.rename(tmpname, self.filename)

class Stopwatch(object):
    """Keeps track of execution time"""
    def __init__(self):
        """Start the stopwatch"""
        self.timestamp_start = datetime.datetime.now()
        self.started = self.timestamp_start.strftime('%c')
        self.laptimes = []
    def lap(self):
        """Mark completion of a lap (or cycle, etc.)"""
        self.laptimes.append(datetime.datetime.now())
    def elapsed(self):
        """Returns datetime.timedelta instance of time since start.
        
        The result looks nice when converted to a string.
        """
        return datetime.datetime.now() - self.timestamp_start
    def mean_lap_time(self):
        """Mean time, in seconds, between laps"""
        if not self.laptimes:
            return np.nan
        return (self.laptimes[-1] - self.timestamp_start).total_seconds() \
                / float(len(self.laptimes))
    def estimate_completion(self, total_laps, laps=None):
        """Estimate how much time is left until completion.

        Specify expected total number of laps.
        
        Returns datetime.timedelta, which can be converted to a string.
        Returns 'None' if time would be negative.
        """
        if laps is None:
            laps = len(self.laptimes)
        if laps > total_laps or laps == 0:
            return None
        else:
            return datetime.timedelta(0, self.elapsed().total_seconds() \
                                         / laps * (total_laps - laps))

def _format_td(timedelt):
    """Format a timedelta object as hh:mm:ss"""
    if timedelt is None:
        return ''
    s = int(round(timedelt.total_seconds()))
    hours = s // 3600
    minutes = (s % 3600) // 60
    seconds = (s % 60)
    return '{:02d}:{:02d}:{:02d}'.format(hours, minutes, seconds)

class Progress(StatusFile):
    def __init__(self, persistent_info=None,
                 filename='taskerstatus.json'):
        """persistent_info : dict of information to always report
        """
        persistent_info = persistent_info.copy()
        persistent_info['pid'] = os.getpid()
        super(Progress, self).__init__(persistent_info, filename)
        self.total_parts = None
        self.stopwatch = Stopwatch()
        self.total = None # total units of work, if we ever find out

    def update(self, newinfo):
        newinfo.update({'started': self.stopwatch.started,
                        'elapsed_time': _format_td(self.stopwatch.elapsed())})
        super(Progress, self).update(newinfo)

    def working(self, current=None, total=None, info=None):
        """Report progress on a task.

        current : current unit of work, counting from 0 (optional)
        total : total units of work (optional)
        info : dict of extra information (optional)
        """
        tmpinfo = {'status': 'working'}
        if current is not None:
            tmpinfo['current'] = current + 1
        if total is not None:
            self.total = total
            tmpinfo['total'] = total
        if (current is not None) and current > 0 \
                and (total is not None) and total > 0:
            tmpinfo['time_left'] = _format_td(
                self.stopwatch.estimate_completion(total, current))
        if info is not None:
            tmpinfo.update(info)
        self.update(tmpinfo)

    def tally(self, iterable, total=None, info=None):
        """Pass-through generator that tracks progress.

        total : total units of work (optional)
        info : dict of extra information to report (optional)

        Example:
        for frame in tally(frames, len(frames)):
            do_something_with(frame)
        """
        for i, item in enumerate(iterable):
            self.stopwatch.lap()
            tmpinfo = {'status': 'working',
                'current': i + 1,
                'time_per': self.stopwatch.mean_lap_time()}
            if total is not None:
                self.total = total
                tmpinfo['total'] = total
                if total > 0:
                    tmpinfo['time_left'] = _format_td(
                        self.stopwatch.estimate_completion(total))
            if info is not None:
                tmpinfo.update(info)
            self.update(tmpinfo)
            yield item

    def finish(self, info=None):
        """Signal end of task.

        info : dict of extra info to report (optional)

        (This is called automatically by Task.)
        """
        tmpinfo = {'status': 'done'}
        if self.total is not None:
            tmpinfo['total'] = self.total
        if info is not None:
            tmpinfo.update(info)
        self.update(tmpinfo)
