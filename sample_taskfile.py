# The package will provide a function to (re)load this file and extract this function.
from tasker import TaskMan, JSON, Pandas
import pandas

class MovieTasker(TaskMan):
    # This could have @cachedprop attributes, basic metadata, etc.
    pass
def use(dirname):
    task = MovieTasker(dirname)
    # Config may be set right here, or by MovieTasker.__init__(), or both.
    task.conf = dict(one='one', two=2.0, name=task.name)
    # We define the tasks here so that we can use closures for maximum
    # directory-specific flexibilty
    @task([], JSON('one.json'))
    def one(ins):
        """Returns string set in config."""
        return task.conf['one'] # Goes to JSON
    @task(one, [JSON('two.json'), JSON('2b.json')])
    def two(input):
        return task.conf['two'], {'twofloat': task.conf['two'], 
                'onestr': input, 'name': task.conf['name']}
    @task([one, two], Pandas('three.h5'))
    def three(ins):
        """Returns a Pandas Series"""
        twofloat = ins[1][1]['twofloat']
        return pandas.Series([twofloat,])
    #... now available as task.one(), task.one.refresh(), etc.
    return task # Very important!
