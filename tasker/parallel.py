import itertools

def imap_throttled_ipyparallel(func, iterable, load_balanced_view=None, 
                               wait_interval=0.05, buffer_factor=3):
    """Asynchronous parallel "imap()" with bounded buffer.
    
    The accumulation of computed results in memory is limited
    to twice the number of nodes.
    
    Parameters:
        func : function
            func(x) will be called for each input item, as in map().
        iterable : iterable
        load_balanced_view : ipyparallel.LoadBalancedView, optional
            The cluster to be used for the computation. If
            omitted or "None", the computation is performed
            locally.
        wait_interval : float
            Time in seconds to wait for new results to become ready.
        buffer_factor : positive integer, optional
            For each node in the cluster, allow at most this many
            completed results in memory.
    
    Yields:
        Results of func() in order of "iterable", as they become
        available.
        
    Notes:
        "wait_interval"/"buffer_factor" should be shortened/increased
        if the execution time of "func" is very heterogeneous.
    """
    # Plan: Cram the queue up to its maximum size.
    # Wait for the earliest-submitted task to finish, then
    # yield its result. As spots become available, add more tasks.
    # Once the input has been exhausted, wait on the remaining tasks.
    if not isinstance(buffer_factor, int) or buffer_factor < 1:
        raise ValueError('buffer_factor must be a positive integer')
        
    if load_balanced_view is None:
        for item in iterable:
            yield func(item)
        raise StopIteration
        
    queue = []  # deque would be marginally faster.
    for item in iterable:
        queue.append(load_balanced_view.apply_async(func, item))
        while len(queue) >= len(load_balanced_view) * buffer_factor:
            queue[0].wait(wait_interval)  # Returns instantly if ready
            if queue[0].ready():
                yield queue.pop(0).get()
    # Wait on the outstanding tasks
    for result in queue:
        yield result.get()
