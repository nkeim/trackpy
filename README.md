trackpy
=======

An impure python implementation of Crocker-Grier for single particle tracking.

Differs from the original (Caswell) version mostly in that it requires the [`numba`](http://numba.pydata.org) Python-to-LLVM dynamic compiler to dramatically speed up some key algorithms.

Installation
----

### Easiest

The easiest route is to first install the free [Anaconda](https://store.continuum.io/cshop/products/) Python distribution from [Continuum Analytics](http://continuum.io). 

You can install directly from this repository, no download necessary, with

    pip install 'git+http://github.com/nkeim/trackpy/@numba#egg=trackpy'

*Tip for novices: Be sure that the `pip` or `python` you are running in these directions belongs to the Python installation you'll use for tracking (e.g. anaconda).*

This also works for upgrading to the latest version. If you want to muck about with the source code, add a `-e` after `pip install`, and a `src` directory will be created in your current directory.

### The other way

Download the source and run

    python setup.py install

in the source directory. 

If you're not using Anaconda, make sure you have the following Python packages:

- `numba`
- `scipy` and `numpy`


Documentation
---

[Documentation](http://tacaswell.github.com/trackpy/)
