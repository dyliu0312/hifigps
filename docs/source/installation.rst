Installation
=============

Requirements
------------

- Python >= 3.8
- numpy
- scipy
- matplotlib
- astropy
- h5py

Optional dependencies:

- ``tqdm`` — for progress bars in ``hifigps-stack``
- ``illustris_python`` — for ``hifigps-find-fuzzy``

Install from source
-------------------

Clone and install::

    git clone https://github.com/dyliu0312/hifigps.git
    cd hifigps
    pip install .

With test dependencies::

    pip install .[test]

With all optional dependencies::

    pip install .[test,prepare,stack]
