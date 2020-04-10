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

import numpy as np
import sparse

import fwdpy11

# Functions related to FS calculation

NOT_A_SAMPLE = np.iinfo(np.int32).min


def _include_both(m):
    return True


def _include_neutral(m):
    return m.neutral is True


def _include_selected(m):
    return m.neutral is False


def _validate_windows(windows, genome_length):
    if windows != sorted(windows, key=lambda x: x[0]):
        raise ValueError("windows must be sorted in increasing order")
    for w in windows:
        if w[0] >= w[1]:
            raise ValueError("windows must be [a, b) and a < b")
        for i in w:
            if i < 0 or i > genome_length:
                raise ValueError("window coordinates must be [0, genome_length)")
    for i in range(1, len(windows), 2):
        if windows[i - 1][0] < windows[i][1] and windows[i][0] < windows[i - 1][1]:
            raise ValueError("windows cannot overlap")


def _update_window(windex, tree, windows):
    while windex < len(windows) and tree.right > windows[windex][1]:
        windex += 1
    return windex


def _position_in_window(pos, window):
    return pos >= window[0] and pos < window[1]


def _simplify(tables, samples, simplify):
    if simplify is False:
        return tables, samples

    t, i = fwdpy11.simplify_tables(tables, samples)
    return t, i[samples]


def _1dfs(self, samples, windows, include_function, simplify):
    """
    Returns an array with the zero and fixed
    bins masked out.  The masking is for
    consistency w/ndfs output.
    """
    t, s = _simplify(self, samples, simplify)
    fs = [np.zeros(len(s) + 1, dtype=np.int32) for i in windows]
    ti = fwdpy11.TreeIterator(t, s)
    windex = 0
    sites = np.array(t.sites, copy=False)
    positions = sites["position"][:]
    for tree in ti:
        for m in tree.mutations():
            if include_function(m):
                pos = positions[m.site]
                while windex < len(windows) and windows[windex][1] < pos:
                    windex += 1
                if windex >= len(windows):
                    break
                if _position_in_window(pos, windows[windex]):
                    c = np.uint32(tree.leaf_counts(m.node))
                    fs[windex][c] += 1

    rv = []
    for i in fs:
        i = np.ma.array(i)
        i[0] = np.ma.masked
        i[-1] = np.ma.masked
        rv.append(i)
    return rv


def _ndfs(
    self, samples, sample_groups, num_sample_groups, windows, include_function, simplify
):
    """
    For more than one sample, always work with the joint FS
    at first. When the marginals are wanted, we extract them
    later.
    """
    shapes = tuple(len(i) + 1 for i in samples)
    dok_JFS = [sparse.DOK(shapes, dtype=np.int32) for i in windows]

    sample_list = np.where(sample_groups != NOT_A_SAMPLE)[0]
    t, s = _simplify(self, sample_list, simplify)
    ti = fwdpy11.TreeIterator(t, s, update_samples=True)
    counts = np.zeros(len(samples), dtype=np.int32)
    windex = 0
    sites = np.array(t.sites, copy=False)
    positions = sites["position"][:]
    for tree in ti:
        for m in tree.mutations():
            if include_function(m):
                pos = positions[m.site]
                while windex < len(windows) and windows[windex][1] < pos:
                    windex += 1
                if windex >= len(windows):
                    break
                if _position_in_window(pos, windows[windex]):
                    counts[:] = 0
                    d = tree.samples_below(m.node)
                    if len(d) > 0:
                        c = np.unique(sample_groups[d], return_counts=True)
                        counts[c[0]] += c[1]
                        dok_JFS[windex][tuple((i) for i in counts)] += 1
    return [sparse.COO(i) for i in dok_JFS]


def _handle_fs_marginalizing(fs, marginalize, nwindows, num_sample_groups):
    if marginalize is False or num_sample_groups == 1:
        return fs
    if nwindows == 1:
        temp = {}
        for i in range(num_sample_groups):
            axes = tuple(j for j in range(num_sample_groups) if j != i)
            temp[i] = np.ma.array(fs.sum(axis=axes).todense())
            temp[i][0] = np.ma.masked
            temp[i][-1] = np.ma.masked
        return temp

    rv = []
    for f in fs:
        temp = {}
        for i in range(num_sample_groups):
            axes = tuple(j for j in range(num_sample_groups) if j != i)
            temp[i] = np.ma.array(f.sum(axis=axes).todense())
            temp[i][0] = np.ma.masked
            temp[i][-1] = np.ma.masked
        rv.append(temp)
    return rv


def _fs_implementation(self, samples, windows, include_function, simplify):
    """
    self is a fwdpy11.TableCollection
    samples is a list of 1d numpy.ndarray
    """
    if len(samples) == 0:
        raise ValueError("empty list of samples")

    # Don't figure out sample groups when
    # there is only 1
    if len(samples) == 1:
        return _1dfs(self, samples[0], windows, include_function, simplify)

    sample_groups = np.array([NOT_A_SAMPLE] * len(self.nodes), dtype=np.int32)
    for i, j in enumerate(samples):
        if np.any(sample_groups[j] != NOT_A_SAMPLE):
            raise ValueError("sample nodes cannot be part of multiple sample groups")
        sample_groups[j] = i

    num_sample_groups = i + 1
    return _ndfs(
        self,
        samples,
        sample_groups,
        num_sample_groups,
        windows,
        include_function,
        simplify,
    )


def _fs(
    self,
    samples,
    marginalize=False,
    simplify=False,
    windows=None,
    separate_windows=False,
    include_neutral=True,
    include_selected=True,
):
    """
    :param samples: lists of numpy arrays of sample nodes
    :param marginalize: For ``FS`` involving multiple samples,
                        extract out the marginal ``FS`` per sample.
    :param simplify: If ``True``, simplify with respect to the sample
                     set prior to calculating the ``FS``.
    :param windows: A list of non-overlapping intervals from which
                    the ``FS`` is calculated.
    :param separate_windows: If ``True``, return ``FS`` separately for
                             each interval in ``windows``.
    :param include_neutral: Include neutral mutations?
    :param include_selected: Include selected mutations?

    :returns: The mutation frequency spectrum
    :rtype: object

    The details of the return value depend heavily on the options.

    * If a single sample list is provided, the return value is a
      1-d ndarray with the first and last bins masked.  The unmasked
      bins correspond to frequencies 1 to n-1, where n is the number
      of sample nodes.
    * If multiple sample lists are provided, the return value is a
      ``sparse.COO`` matrix.  For each dimension, the 0 and n bins
      are included.
    * If ``marginalize == True`` and more than one sample bin is provided,
      the ``sparse.COO`` matrix is converted into a dict where the key
      is the index of the sample list and the value is a dense 1-d array,
      masked as described above.
    * If multiple windows are provided and ``separate_windows == False``,
      then the return value is a single frequency spectrum summed over windows.
      If ``separate_windows == True``, then a list of frequency spectra is
      returned, indexed in the same order as the input windows.

    .. versionadded:: 0.6.0

        Python implementation added
    """
    for s in samples:
        if len(s) == 0:
            raise ValueError("samples lists cannot be empty")
        if len(s) < 2:
            raise ValueError("samples lists must have at least two nodes")
        if np.any(s >= len(self.nodes)):
            raise ValueError("invalid samples")

    if windows is None:
        windows = [(0, self.genome_length)]

    _validate_windows(windows, self.genome_length)

    if include_neutral is False and include_selected is False:
        raise ValueError(
            "One or both of include_neutral " "and include_selected must be True"
        )

    include_function = _include_both
    if include_neutral is False:
        include_function = _include_selected
    elif include_selected is False:
        include_function = _include_neutral

    lfs = _fs_implementation(self, samples, windows, include_function, simplify)
    if separate_windows is True:
        lfs = _handle_fs_marginalizing(lfs, marginalize, len(windows), len(samples))
        return lfs
    fs = lfs[0]
    for i in lfs[1:]:
        fs += i
    fs = _handle_fs_marginalizing(fs, marginalize, len(windows), len(samples))
    return fs


def _patch_table_collection(t):
    t.fs = _fs
