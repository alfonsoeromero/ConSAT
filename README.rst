======
ConSAT
======

.. image:: https://travis-ci.org/paccanarolab/ConSAT.svg?branch=hackaton
    :target: https://travis-ci.org/paccanarolab/ConSAT

.. image:: https://api.codeclimate.com/v1/badges/46c5df805d70f2dec2f8/maintainability
   :target: https://codeclimate.com/github/paccanarolab/ConSAT/maintainability
   :alt: Maintainability

.. image:: https://api.codeclimate.com/v1/badges/46c5df805d70f2dec2f8/test_coverage
   :target: https://codeclimate.com/github/paccanarolab/ConSAT/test_coverage
   :alt: Test Coverage

-----------------------------------------
The Consensus Signature Architecture Tool
-----------------------------------------

:Author: Alfonso E. Romero, Tamás Nepusz, Rajkumar Sasidharan, Alberto Paccanaro

This is the documentation of ``consat``, a Python-based application to provide protein
families defined as consensus domain architectures. ``consat`` is the stand-alone tool
used in the homonym web server (`http://www.paccanarolab.org/consat`_).

.. _`http://www.paccanarolab.org/consat`: http://www.paccanarolab.org/consat

This tool takes a file of protein sequences, a file with domain assignments in InterProScan
format, and a file with a library of HMMs domains not included in InterPro, and produces
an assignment of proteins to protein families, as well as functional information association
(mainly GO terms) based in those families.

As this software is an evolution of ``gfam`` (and also includes it), most of the changes
have been propagated backwards, and this gives the possibility of running ``gfam`` with the
new developments implemented in ``consat`` as well as with improvements in the original 
algorithm.

Requirements
============

You will need the following tools to run ``consat``:

* `Python 3.6`_. 

* `HMMer 3.0`_ (or newer versions).

* `NCBI BLAST`_; if you are using ``gfam`` from this software package. 

.. _`Python 3.6`: http://www.python.org
.. _`HMMer 3.0`: http://hmmer.janelia.org
.. _`NCBI BLAST`: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST

The latest release of `SciPy`_ is recommended, but not necessary.
``gfam`` uses `SciPy`_ for calculating the logarithm of the gamma
function in the overrepresentation analysis routines, but it falls
back to a (somewhat slower) Python implementation if `SciPy`_ is
not installed.

.. _`SciPy`: http://www.scipy.org

Running ``consat``
==================

``consat`` is driven by a master configuration file named ``consat.cfg``.
A sample configuration file is given in the distribution. The sample
file works fine for the gene sequences of *Arabidopsis thaliana*; for
other species, you might have to tweak some of the parameters, and you
will surely have to modify the paths to the data files and maybe other
parameters. The configuration file is documented and mostly 
self-explanatory.

You can launch ``consat`` by typing::

    $ bin/consat

This will run the whole ``consat`` analysis pipeline using the configuration
specified in ``consat.cfg``. If your configuration file is named otherwise,
you can run it by typing::

    $ bin/consat -c my_config.cfg

Questions, comments
===================

If you have a question or a comment about ``consat`` or you think you have
found a bug, feel free to `contact me`_.

.. _contact me: http://www.cs.rhul.ac.uk/~aeromero
