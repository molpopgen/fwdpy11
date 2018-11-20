.. _binary_pops:

Writing populations to files in binary format
====================================================================================

.. versionadded:: 0.2.0

Beginning with version 0.2.0, fwdpy11 allows population objects to be written directly 
to binary-format files and read back in from such files.  This funcitonality is carried
out by the following functions:

* :func:`fwdpy11.SlocusPop.dump_to_file`
* :func:`fwdpy11.MlocusPop.dump_to_file`
* :func:`fwdpy11.SlocusPop.load_from_file`
* :func:`fwdpy11.MlocusPop.load_from_file`

The writing functions are limited to write-only, no compression, and no support for appending to a file.

