### Everything below will be in a project-specific file, taskfile.py.
# The package will provide a function to (re)load this file and extract this function.

class MovieTasker(TaskMan):
    # This could have @cachedprop attributes, basic metadata, etc.
    # As well as knowing some fun things to do with the 'conf' dict, like 
    # write to disk
    pass
def movie(dirname):
    task = MovieTasker(dirname)
    task.conf = dict(one='one', two=2.0, name=task.name)
    # We define the tasks here so that we can use closures for maximum
    # directory-specific flexibilty
    @task([], 'one.json')
    def one(ins):
        """ONE"""
        print task.conf['one']
        return task.conf['one']
    # As we go through this, we can set additional config items, 
    # conditionally define tasks, etc.
    @task(one)
    def two(ins):
        print task.conf['two']
        print ins[0] # or ins['one:one.json'] or ins['one'], etc.
        return 1
    #... now available as task.one(), task.one.refresh(), etc.
    return task

t = movie('.')
