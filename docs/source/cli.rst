Command Line Interface
=======================

hifigps-convolve
----------------

Convolve a map with the FAST main beam.

.. code-block:: bash

    hifigps-convolve MAP_FILE OUT_PATH NWORKER [-z REDSHIFT] [-k KEY] [-d {float32,float64}]

**Arguments:**

- ``MAP_FILE`` — Path to input map file (.hdf5)
- ``OUT_PATH`` — Path to output directory
- ``NWORKER`` — Number of worker processes
- ``-z, --redshift`` — Redshift (default: 0.1)
- ``-k, --key`` — HDF5 key for map data (default: T)
- ``-d, --dtype`` — Output data type (default: float32)

**Example:**

.. code-block:: bash

    hifigps-convolve /data/test.hdf5 /output/ 4 -z 0.1

hifigps-find-fuzzy
------------------

Find inner-fuzzy particles in IllustrisTNG simulations (see ![TNG docs](https://www.tng-project.org/data/docs/specifications/) for more detials about fuzzy particles).

.. code-block:: bash

    hifigps-find-fuzzy BASE SNAP [-o OUTPUT] [-p {gas,dm,star,bhs,bhw,stars}]

**Arguments:**

- ``BASE`` — Base path to TNG simulation data
- ``SNAP`` — Snapshot number
- ``-o, --output`` — Output file path (required)
- ``-p, --part-type`` — Particle type (default: gas)

**Example:**

.. code-block:: bash

    hifigps-find-fuzzy /path/to/TNG/ 91 -o /output/fuzz.h5

hifigps-stack
--------------

Galaxy pairwise stacking.

.. code-block:: bash

    hifigps-stack --map-base MAP_BASE --map-prefix MAP_PREFIX
                   --paircat-base PAIRCAT_BASE --paircat-prefix PAIRCAT_PREFIX
                   --out-base OUT_BASE --out-prefix OUT_PREFIX
                   --nfs NFS --ssize SSIZE [options]

**Required arguments:**

- ``--map-base`` — Base path to input map
- ``--map-prefix`` — Prefix of input map file
- ``--paircat-base`` — Base path to pair catalog
- ``--paircat-prefix`` — Prefix of pair catalog file
- ``--out-base`` — Base path to output
- ``--out-prefix`` — Prefix of output file
- ``--nfs`` — Number of frequency slices to stack
- ``--ssize`` — Split size for processing

**Optional arguments:**

- ``--nworker`` — Number of workers (default: CPU count)
- ``--random-flip`` — Randomly flip signal (default: True)
- ``--halfwidth`` — Stack result map half-width (default: 3.0)
- ``--npix-x``, ``--npix-y`` — Stack result map pixels (default: 120)
- ``--compression`` — gzip/lzf/none (default: gzip)
- ``--skip-exist`` — Skip existing splits (default: False)

**Example:**

.. code-block:: bash

    hifigps-stack --map-base /data/ --map-prefix map \
                  --paircat-base /data/ --paircat-prefix paircat \
                  --out-base /output/ --out-prefix stack \
                  --nfs 10 --ssize 1000 --nworker 24
