Orange
======

[![build: passing](https://img.shields.io/travis/biolab/orange3.svg)](https://travis-ci.org/biolab/orange3)
[![codecov](https://codecov.io/gh/biolab/orange3/branch/master/graph/badge.svg)](https://codecov.io/gh/biolab/orange3)

[Orange] is a component-based data mining software. It includes a range of data
visualization, exploration, preprocessing and modeling techniques. It can be
used through a nice and intuitive user interface or, for more advanced users,
as a module for the Python programming language.

This is a development version of Orange 3. The stable version 2.7 is still
available ([binaries] and [sources]).

[Orange]: http://orange.biolab.si/
[binaries]: http://orange.biolab.si/orange2/
[sources]: https://github.com/biolab/orange


Installing
----------
This version of Orange requires Python 3.4 or newer. To build it and install
it in a development environment, run:

    # Install some build requirements via your system's package manager
    sudo apt-get install virtualenv git build-essential

    # Also install Qt dependencies for the GUI
    sudo apt-get install python3-pyqt4
    # or if python version is >= 3.5
    #  pip install pyqt5 

    # Create a separate Python environment for Orange and its dependencies ...
    virtualenv --python=python3 --system-site-packages orange3venv
    # ... and make it the active one
    source orange3venv/bin/activate

    # Clone the repository and move into it
    git clone https://github.com/biolab/orange3.git
    cd orange3

    # Install the minimum required dependencies first
    pip install -r requirements-core.txt  # For Orange Python library
    pip install -r requirements-gui.txt   # For Orange GUI

    pip install -r requirements-sql.txt   # To use SQL support
    pip install -r requirements-opt.txt   # Optional dependencies, may fail

    # Finally install Orange in editable/development mode.
    pip install -e .

Installation of SciPy and qt-graph-helpers is sometimes challenging because of
their non-python dependencies that have to be installed manually. More
detailed, if mostly obsolete, guides for some platforms can be found in
the [wiki].

[wiki]: https://github.com/biolab/orange3/wiki

Anaconda Installation
---------------------

First, install [Anaconda] for your OS (Python version 3.5+). Create virtual environment for Orange:

    conda create python=3 --name orange3 

In your Anaconda Prompt add conda-forge to your channels:

    conda config --add channels conda-forge

This will enable access to the latest Orange release. Then install Orange3:

    conda install orange3

[Anaconda]: https://www.continuum.io/downloads

Starting Orange GUI
-------------------

Orange GUI requires PyQt, which is not pip-installable in Python 3. You
have to download and install it system-wide. Make sure that the virtual
environment for orange is created with `--system-site-packages`, so it will
have access to the installed PyQt4.

To start Orange GUI from the command line, assuming it was successfully
installed, run:

    orange-canvas
    # or
    python3 -m Orange.canvas

Append `--help` for a list of program options.

If you're running Orange with PyQt5 or if you have multiple PyQt versions
available, set the environmental variable `QT_API` to the PyQt version to use,
e.g.:

    export QT_API=pyqt5
    orange-canvas


Windows dev setup
-----------------

Windows + GCC:

    python setup.py build_ext --inplace --compile=mingw32

Evolution of genes
-----------------

### New features
#### Multiple Sequence Alignment

The Multiple Sequence Alignment widget takes in multiple sequences and calculates minimal distance and respective alignment for each pair of sequences.
It implements Needleman-Wunsch dynamic programming algorithm for optimal global alignment.

Widget enables setting custom scores:
* gap score - a non negative integer *(default = 1)*
* substitution score - a non negative integer *(default = 1)*
* alignment score - zero or negative integer *(default = 0)*


	inputs = [("Data", Orange.data.Table, "set_data")]
    outputs = [("Distances", Orange.misc.DistMatrix),
               ("Alignments", Orange.data.Table)]


#### Alignment Distance Matrix

Alignment Distance Matrix widget takes in pairwise distances and alignments from Multiple Sequence alignment widget and
displays them in a manner, similar to that of the Distance Matrix widget. Clicking on the cell corresponding to a pair
of sequences displays their optimal alignment.

	inputs = [("Distances", DistMatrix, "set_distances"),
              ("Alignments", Table, "set_alignments")]
    outputs = [("Distances", DistMatrix),
               ("Table", Table)]
### Case study

The image below shows a minimal pipeline for obtaining and visualizing sequence alignments using Multiple Sequence
Alignment widget.

1. Data is input using the *File* widget. Sequence names are *meta* type strings and sequences' type has
to be *nominal* and is set as *feature* in the widget's dialog window.

2. Rows of interest are selected using the *Data Table* widget.

3. Selected data is sent to the *Multiple Sequence Alignment* widget. Score values can be set using the widget's dialog
or can be left at default values.

4. Distances and alignments are passed on to the *Alignment Distance Matrix* widget. Its dialog shows a table of all distances.
Sequence names can be displayed for improved readability. Clicking on a cell in the table displays optimal alignment for
the chosen pair of sequences at the bottom of the dialog.


 ![enter image description here](https://lh3.googleusercontent.com/-I3_khTVbleU/WHHtVOflGyI/AAAAAAAAADI/MUDmVVULPOUNX9eGQ1uKPosFdpNbU_99wCLcB/s0/workflow.PNG "workflow.PNG")

