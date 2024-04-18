from typing import Iterable

import numpy as np

from .._fwdpy11 import Edge, MutationRecord, Node, Site, ll_TableCollection

# Functions related to FS calculation

NOT_A_SAMPLE = np.iinfo(np.int32).min


class TableCollection(ll_TableCollection):
    def __init__(self, *args):
        super(TableCollection, self).__init__(*args)

    @property
    def L(self) -> float:
        return self._genome_length

    @property
    def genome_length(self) -> float:
        """Genome length"""
        return self._genome_length

    @property
    def edges(self) -> Iterable[Edge]:
        """Supports buffer protocol"""
        return self._edges

    @property
    def nodes(self) -> Iterable[Node]:
        """Supports buffer protocol"""
        return self._nodes

    @property
    def sites(self) -> Iterable[Site]:
        """Supports buffer protocol"""
        return self._sites

    @property
    def mutations(self) -> Iterable[MutationRecord]:
        """Supports buffer protocol"""
        return self._mutations

    @property
    def input_left(self):
        """Edge input order"""
        return self._input_left

    @property
    def output_right(self):
        """Edge output order"""
        return self._output_right

    def build_indexes(self):
        """
        Build edge input/output indexes.

        .. versionadded:: 0.13.0

        """
        self._build_indexes()

    def preserved_nodes(self):
        raise AttributeError(
            "TableCollection.preserved_nodes was removed in 0.13.0 and is"
            " now an attribute of DiploidPopulation"
        )
