from typing import IO, Dict, Iterable, Iterator, List, Optional, Tuple, Union

import demes
import fwdpy11._types
import fwdpy11.tskit_tools._dump_tables_to_tskit
import numpy as np
import tskit  # type: ignore

from .._fwdpy11 import DiploidGenotype, DiploidMetadata, ll_DiploidPopulation
from .model_params import ModelParams
from .new_mutation_data import NewMutationData
from .population_mixin import PopulationMixin
from .table_collection import TableCollection


class DiploidPopulation(ll_DiploidPopulation, PopulationMixin):
    """
    A population of diploid individuals.

    To initialize a population of size N and genome length L:

    fwdpy11.DiploidPopulation(N, L)

    To initialize a population with two demes of initial
    sizes N1 and N2 and genome length L:

    fwdpy11.DiploidPopulation([N1, N2], L)

    .. versionchanged:: 0.16.0

        Added __copy__ and __deepcopy__.
        In general, the correct way to copy
        instances if this class is via :func:`copy.deepcopy`.
        Calling :func:`copy.copy` will copy some of the underlying
        C++ objects, but the two objects will share the same
        table collection.  This situation of sharing is almost certainly
        not what one wants.
    """

    def __init__(
        self,
        N: Union[int, List[int]],
        length: float,
        *,
        ll_pop: Optional[ll_DiploidPopulation] = None,
    ):
        if ll_pop is None:
            super(DiploidPopulation, self).__init__(N, length)
        else:
            super(DiploidPopulation, self).__init__(ll_pop)

        self._pytables = TableCollection(self._tables)

    def __copy__(self):
        ll_pop = super(DiploidPopulation, self).__copy__()
        return self.__class__(0, 0.0, ll_pop=ll_pop)

    def __deepcopy__(self, memo):
        ll_pop = super(DiploidPopulation, self).__deepcopy__(memo)
        return self.__class__(0, 0.0, ll_pop=ll_pop)

    @classmethod
    def create_from_tskit(cls, ts: tskit.TreeSequence,
                          *, import_mutations=False):
        """
        Create a new object from an tskit.TreeSequence

        :param ts: A tree sequence from tskit
        :type ts: tskit.TreeSequence
        :param import_mutations: if `True`, import mutations and sites
                                 from the tree sequence
        :type import_mutations: bool

        :return: A population object with an initialized
                 :class:`fwdpy11.TableCollection`
        :rtype: :class:`fwdpy11.DiploidPopulation`

        .. note::

            In general, initializing a population using
            the output from a coalescent simulation is
            a tricky business.  There are issues of
            parameter scaling and the appropriateness
            of the coalescent model itself. A key issue
            is that your input tree sequence must have
            node times in the correct time units! (Generations,
            for example.)

        .. warning::

            Support for importing mutations is currently limited to
            two cases:

            1. The tree seq has mutations entirely from tskit's API
               (see :ref:`import_mutations_from_tskit_vignette`).
            2. The tree seq has mutations entirely from a previous run
               of fwdpy11.


        .. versionadded:: 0.2.0

        .. versionchanged:: 0.19.8

           Correctly handle two cases:

           1. The tree seq has mutations entirely from tskit's API
              (see :ref:`import_mutations_from_tskit_vignette`).
           2. The tree seq has mutations entirely from a previous run
              of fwdpy11.

        """
        generation = 0
        if len(ts.metadata) > 0:
            if "generation" in ts.metadata:
                generation = ts.metadata["generation"]
        ll = ll_DiploidPopulation._create_from_tskit(ts, generation)
        if import_mutations is True:
            import fwdpy11.tskit_tools
            for mutation in ts.mutations():
                # if mutation.metadata["origin"] < 0.0:
                #     msg = "all origin fields in mutation metadata"
                #     msg += "must be non-negative"
                #     raise ValueError(msg)
                if not mutation.time.is_integer():
                    raise ValueError("mutation times cannot contain decimals")
            mutation_nodes = [i.node for i in ts.mutations()]
            decoded_md = fwdpy11.tskit_tools.decode_mutation_metadata(ts)
            try:
                # The case of the tree seq has all mutations added
                # to it via tskit API
                for i, j in zip(ts.mutations(), decoded_md):
                    if j is not None:
                        if i.time != j.g:
                            msg = f"metadata origin time {j.g} != "
                            msg += f"mutation time {i.time} in tables"
                            raise ValueError(msg)
            except ValueError:
                # The case of the tree seq has all of its mutations
                # from a pervious fwdpy11 run
                keys = set([i.key for i in decoded_md])
                for i in ts.mutations():
                    md = i.metadata
                    key = (ts.site(i.site).position, md['s'], md['origin'])
                    if key not in keys:
                        raise ValueError(f"{key}")

            mvec = fwdpy11._fwdpy11.MutationVector()
            for i in decoded_md:
                if i is not None:
                    mvec.append(i)
            ll._set_mutations(mvec, mutation_nodes,
                              [int(i.time) - ll._generation for i in ts.mutations()])
        return cls(0, 0.0, ll_pop=ll)

    @ classmethod
    def load_from_file(cls, filename: str):
        """
        Load population from the output of
        :func:`fwdpy11.DiploidPopulation.dump_to_file`.

        :param filename: A file name
        :type filename: str
        """
        ll = ll_DiploidPopulation._load_from_file(filename)
        return cls(0, 0.0, ll_pop=ll)

    @ classmethod
    def load_from_pickle_file(cls, filename: IO):
        """
        Read in a pickled population from a file.
        The file muse have been generated by
        a call to :func:`fwdpy11.DiploidPopulation.pickle_to_file`.

        :param f: A handle to a file opened in 'rb' mode.

        .. versionadded: 0.3.0

        """
        ll = ll_DiploidPopulation._load_from_pickle_file(filename)
        return cls(0, 0.0, ll_pop=ll)

    @ property
    def alive_nodes(self) -> np.ndarray:
        """
        List of alive nodes corresponding to individuals.

        .. versionadded:: 0.5.3
        """
        md = np.array(self.diploid_metadata, copy=False)
        return md["nodes"].flatten()

    @ property
    def ancient_sample_nodes(self) -> np.ndarray:
        """
        Return an array of nodes associated with preserved/ancient samples.

        Alias for :attr:`fwdpy11.DiploidPopulation.preserved_nodes`.

        .. versionadded:: 0.13.0
        """
        return self.preserved_nodes

    @ property
    def ancient_sample_metadata(self) -> Iterable[DiploidMetadata]:
        """Supports buffer protocol"""
        return self._ancient_sample_metadata

    @ property
    def diploids(self) -> Iterable[DiploidGenotype]:
        """Supports buffer protocol"""
        return self._diploids

    @ property
    def diploid_metadata(self) -> Iterable[DiploidMetadata]:
        """Supports buffer protocol"""
        return self._diploid_metadata

    @ property
    def preserved_nodes(self) -> np.ndarray:
        """
        Return an array of nodes associated with preserved/ancient samples.

        :rtype: numpy.ndarray

        .. versionadded:: 0.13.0

        """
        return np.array(self.ancient_sample_metadata, copy=False)["nodes"].flatten()

    @ property
    def tables(self) -> TableCollection:
        """Access the :class:`fwdpy11.TableCollection`"""
        return self._pytables

    def _get_times(self):
        amd = np.array(self.ancient_sample_metadata, copy=False)
        nodes = np.array(self.tables.nodes, copy=False)
        times = nodes["time"][amd["nodes"][:, 0]]
        utimes = np.unique(times)
        return times, utimes

    def deme_sizes(self, as_dict=False):
        """
        Return the number of individuals in each deme.

        :param as_dict: If True, return results as a dict.
                        If False, numpy arrays are returned

        The default behavior of this function is equivalent to:

        .. code-block:: python

            md = numpy.array(self.diploid_metadata, copy=False)
            return numpy.unique(md['deme'], return_counts=True)

        .. versionadded:: 0.6.0
        """
        md = np.array(self.diploid_metadata, copy=False)
        deme_sizes = np.unique(md["deme"], return_counts=True)
        if as_dict is False:
            return deme_sizes
        return {i: j for i, j in zip(deme_sizes[0], deme_sizes[1])}

    def dump_tables_to_tskit(
        self,
        *,
        model_params: Optional[Union[ModelParams,
                                     Dict[str, ModelParams]]] = None,
        demes_graph: Optional[demes.Graph] = None,
        population_metadata: Optional[Dict[int, object]] = None,
        data: Optional[object] = None,
        seed: Optional[int] = None,
        parameters: Optional[Dict] = None,
        destructive=False,
    ) -> tskit.TreeSequence:
        """
        Dump the population's TableCollection into
        an tskit TreeSequence

        :param model_params: Model parameters to be stored as top-level metadata
        :type model_params: :class:`fwdpy11.ModelParams` or :class:`dict`

        :param demes_graph: A demographic model specified via `demes`.
        :type demes_graph: :class:`demes.Graph`

        :param population_metadata: A mapping from integer id of a deme/population to metadata
        :type population_metadata: dict

        :param data: User-generated data to add to top-level metadata
        :type data: object

        :param seed: Random number seed for top-level metadata.
                     If not `None`, `seed` must be `>= 0`.
        :type seed: int

        :param parameters: The simulation parameters for the provenance table.
        :type parameters: None or dict

        :param destructive: If `True`, delete data held by the current instance.
        :type destructive: bool

        :rtype: tskit.TreeSequence

        For examples, see :ref:`tskitconvert_vignette`.

        .. warning::

            If `destructive` is `True`, further opertations on the
            current instance should be considered undefined behavior
            that may lead to a crash.

        .. versionchanged:: 0.8.2

            Added `parameters`.
            Generate provenance information for return value.
            The provenance information is validated using
            :func:`tskit.validate_provenance`, which may
            raise an exception.

        .. versionchanged:: 0.10.0

            Use tskit metadata schema.
            Mutation time is now stored in the tskit.MutationTable column.
            Origin time of mutations is part of the metadata.

        .. versionchanged:: 0.14.0

            Added `destructive` option.

        .. versionchanged:: 0.15.0

            Added `model_params`, `demes_graph`, `population_metadata`,
            `data` keyword args.

        .. versionchanged:: 0.15.2

            Fixed bug that could generate a :class:`tskit.PopulationTable`
            with an incorrect number of rows.

        .. versionchanged:: 0.17.0

            The `wrapped` keyword argument is deprecated.

        .. versionadded:: 0.19.0

            Removed deprecated `wrapped` keyword argument.

        """
        return fwdpy11.tskit_tools._dump_tables_to_tskit._dump_tables_to_tskit(
            self,
            model_params=model_params,
            demes_graph=demes_graph,
            population_metadata=population_metadata,
            data=data,
            seed=seed,
            parameters=parameters,
            destructive=destructive,
        )

    def dump_to_file(self, filename: str):
        """
        Write a population to a file in binary format.
        """
        self._dump_to_file(filename)

    def pickle_to_file(self, filename: IO):
        """
        Pickle the population to an open file.

        This function may be preferred over
        the direct pickling method because it uses less
        memory.  It is, however, slower.

        To read the population back in, you must call
        :func:`fwdpy11.DiploidPopulation.load_from_pickle_file`.

        :param f: A handle to an open file

        .. versionadded:: 0.3.0

        """
        self._pickle_to_file(filename)

    def sample_timepoints(
        self, include_alive=True
    ):
        """
        Return an iterator over all sample time points.
        The iterator yields time, nodes, and metadata.

        :param include_alive: If True, include currently-alive individuals.

        .. versionadded :: 0.5.1

        .. versionchanged:: 0.5.3

           Monkey-patched into pybind11 class
        """
        amd = np.array(self.ancient_sample_metadata, copy=False)
        amd.flags.writeable = False
        times, utimes = self._get_times()
        for uti in utimes:
            mdidx = np.where(times == uti)[0]
            sample_nodes = amd["nodes"][mdidx].flatten()
            sample_nodes.flags.writeable = False
            mdslice = amd[mdidx][:]
            mdslice.flags.writeable = False
            yield uti, sample_nodes, mdslice

        if include_alive is True:
            md = np.array(self.diploid_metadata, copy=False)
            nodes = md["nodes"].flatten()
            md.flags.writeable = False
            nodes.flags.writeable = False
            yield self.generation, nodes, md

    def add_mutation(
        self,
        rng: fwdpy11.GSLrng,
        *,
        window: Tuple[float, float] = None,
        ndescendants: int = 1,
        deme: Optional[int] = None,
        data: NewMutationData = None,
    ) -> Optional[int]:
        """
        Add a new mutation to the population.

        .. versionadded:: 0.16.0

        :param rng: Random number generator
        :type rng: fwdpy11.GSLrng

        The following arguments are keyword-only:

        :param window: A window [left, right) within which to place the mutation.
                       The default is `None`, meaning the window is the entire
                       genome.
        :type window: tuple[float, float]
        :param ndescendants: The number of alive nodes carrying the new mutation.
                             Default is `1`, implying that a singleton mutation
                             will be generated.
        :type ndescendants: int
        :param deme: The deme in which to place the new mutation
                     The default is `None`, meaning that alive node demes are
                     not considered.
        :type deme: int
        :type data: The details of the new mutation
        :type data: fwdpy11.NewMutationData

        :returns: The key of the new mutation.
                  (The index of the variant in :attr:`DiploidPopulation.mutations`.)

        Implementation details:

        * A set of nodes are found within `window` that are ancestral to
          exactly `ndescendants` alive nodes.
        * From this list of candidate nodes, one is chosen randomly
          to be the node where we place the new mutation.
        * If the node's parent's time is NULL, then the mutations' origin
          time is the first 64 bit integer ancestral to the node time,
          which could be the node time itself.
          Otherwise, a 64 bit integer is chosen uniformly between the
          node time and the node's parent's time.
        * If `deme` is None, any set of `ndescendants` will be considered.
        * If `deme` is :math:`\geq 0`, all `ndescendants` alive nodes must be
          from that deme.
        * If the `deme`/`ndescendants` requirements cannot be satisified,
          the function returns `None`.
        * The genetic values of individuals are not updated by this function.
          Such updates will happen when the population begins evolution
          forwards in time.

        Exceptions:

        :raises ValueError: for bad input.
        :raises RuntimeError: if this function is applied during a simulation
        :raises RuntimeError: if internal checks fail.
        """
        if window is None:
            _window = (0.0, self.tables.genome_length)
        else:
            _window = window

        if deme is None:
            _deme = -1
        else:
            _deme = deme

        if ndescendants is None:
            raise ValueError(f"ndescendants must be > 0: got {ndescendants}")

        from fwdpy11._fwdpy11 import _add_mutation

        key = _add_mutation(
            rng, _window[0], _window[1], ndescendants, _deme, data, self
        )
        if key == np.iinfo(np.uint64).max:
            return None
        return key
