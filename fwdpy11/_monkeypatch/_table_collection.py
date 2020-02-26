#
# Copyright (C) 2020 Kevin Thornton <krthornt@uci.edu>
#
# This file is part of fwdpy11.
#
# fwdpy11 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# fwdpy11 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with fwdpy11.  If not, see <http://www.gnu.org/licenses/>.
#

import fwdpy11
import sparse
import scipy.sparse
import numpy as np

# Functions related to FS calculation


def _1dfs(self, samples):
    """
    Returns an array with the zero and fixed
    bins masked out.  The masking is for
    consistency w/ndfs output.
    """
    fs = np.ma.zeros(len(samples)+1, dtype=np.int32)
    ti = fwdpy11.TreeIterator(self, samples)
    for tree in ti:
        for m in tree.mutations():
            c = tree.leaf_counts(m.node)
            fs[c] += 1

    fs[0] = np.ma.masked
    fs[-1] = np.ma.masked

    return fs


def _ndfs(self, samples, sample_groups, num_sample_groups):
    """
    For more than one sample, always work with the joint FS
    at first. When the marginals are wanted, we extract them
    later.
    """
    shapes = tuple(len(i)+1 for i in samples)
    if num_sample_groups == 2:
        dok_JFS = scipy.sparse.dok_matrix(shapes, dtype=np.int32)
        coo_JFS_type = scipy.sparse.coo_matrix
    else:
        dok_JFS = sparse.DOK(shapes, dtype=np.int32)
        coo_JFS_type = sparse.COO

    sample_list = np.where(sample_groups != -1)[0]
    ti = fwdpy11.TreeIterator(self, sample_list, update_samples=True)
    counts = np.zeros(len(samples), dtype=np.int32)
    for tree in ti:
        for m in tree.mutations():
            counts[:] = 0
            d = tree.samples_below(m.node)
            if len(d) > 0:
                for i in d:
                    counts[sample_groups[i]] += 1
                dok_JFS[tuple((i) for i in counts)] += 1
    return coo_JFS_type(dok_JFS)


def _fs_implementation(self, samples):
    """
    self is a fwdpy11.TableCollection
    samples is a list of 1d numpy.ndarray
    """
    if len(samples) == 0:
        raise ValueError("empty list of samples")

    # Don't figure out sample groups when
    # there is only 1
    if len(samples) == 1:
        return _1dfs(self, samples[0])

    sample_groups = np.array([-1]*len(self.nodes), dtype=np.int32)
    for i, j in enumerate(samples):
        if np.any(sample_groups[j] != -1):
            raise ValueError(
                "sample nodes cannot be part of multiple sample groups")
        sample_groups[j] = i

    num_sample_groups = i + 1
    return _ndfs(self, samples, sample_groups, num_sample_groups)


def _fs(self, samples=None, sample_sizes=None,
        by_deme=True, marginalize=False,
        simplify=False, windows=None,
        separate_windows=False,
        include_neutral=True,
        include_selected=True):
    """
    :param samples: lists of numpy arrays of sample nodes
    :param sample_sizes: Number of nodes of each sample
    :param by_deme: If ``samples and ``sample_sizes`` are both
                    ``None``, obtain ``FS`` separately
                    for each entire deme.
    :param marginalize: For ``FS`` involving multiple samples,
                        extract out the marginal ``FS`` per sample.
    :param simplify: If ``True``, simplify with respect to the sample
                     set prior to calculating the ``FS``.
    :param windows: A list of non-overlapping intervals from which
                    the ``FS`` is calculated.
    :param separate_windows: If ``True``, return ``FS`` separately for
                             each interval in ``windows``.
    :param include_neutral: Include neutral mutations
    :param include_selected: Include neutral mutations
    """
    return _fs_implementation(self, samples)


def _patch_table_collection(t):
    t.fs = _fs
