# Tasker

Tasker is a lightweight Python framework for organizing and automating scientific data analysis. If you have lengthy computations involving several steps, over data in many directories, Tasker will help you

- Store your results to disk.
- Build up a complex analysis out of many intermediate steps, saving the results of earlier steps so you don't have to redo them.
- Only re-compute when necessary — when you change an input file or delete an intermediate result, only the computations that depended on it will be redone.
- Write code that works on data in many similar directories.
- Run said code in parallel.

Check out the quick (yet rather comprehensive) demonstration notebook included with the software.

Unit tests may be run with `py.test`.

# Installation

```
pip install http://github.com/nkeim/tasker/zipball/master
```

If you've already installed Tasker and wish to upgrade:

```
pip install --upgrade http://github.com/nkeim/tasker/zipball/master
```

# About

Written by [Nathan Keim](http://www.seas.upenn.edu/~nkeim/), while at the [Penn Complex Fluids Lab](http://arratia.seas.upenn.edu).

Tasker is named after Tasker Street in Philadelphia, Pennsylvania.
