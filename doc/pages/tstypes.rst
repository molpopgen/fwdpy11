.. _ts_data_types:

Data structures related to tree sequences
======================================================================

To start out, let us consider the follwing tree:

.. figure:: ../images/tree.png

        Fig 1: A tree with seven nodes.

This tree is the "marginal" history of a genomic segment covering the half-open interval :math:`[l, r)`.

Following Kelleher *et al.* (2016, PLoS Computational Biology), we can represent the above tree using two tables:

.. figure:: ../images/tables.png

       Fig 2: The node and edges tables corresponding to Fig 1.

We learn two things from Fig 2:

1. Node tables track the birth times of nodes.  Here, we measure time as increasing from past to the present.
2. Edge tables record the transmissions of genomic intervals from parents to children.  The parent/child fields
   are indexes of the node table.

Edge tables have specific sorting requirements.  The sorting is nested:

1. Decreasing order of parent birth times (as we read the table from top to bottom).
2. For edges with the same parent, child indexes are sorted in increasing order
3. Finally, edges are sorted by increasing left position.

Low-level data types
----------------------------------------------------------

* :class:`fwdpy11.ts.Node` defines nodes
* :class:`fwdpy11.ts.Edge` defines edges
* :class:`fwdpy11.ts.MutationRecord` defines mutations locations on trees and in :class:`fwdpy11.Population` objects.

* :class:`fwdpy11.ts.NodeTable` represents a node table
* :class:`fwdpy11.ts.EdgeTable` represents an edge table
* :class:`fwdpy11.ts.MutationTable` represents a mutation table

Table collections
----------------------------------------------------------

The above data types are encapsulated into the Python class :class:`fwdpy11.ts.TableCollection`.  Instances of this
class are data fields of populations, via :attr:`fwdpy11.Population.tables`.

