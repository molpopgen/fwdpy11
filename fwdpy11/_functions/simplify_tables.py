from typing import List, Tuple, Union

import fwdpy11._fwdpy11
import fwdpy11._types
import numpy as np


def simplify(pop, samples):
    """
    Simplify a TableCollection stored in a Population.

    :param pop: A :class:`fwdpy11.PopulationBase`
    :param samples: A list of samples (node indexes).

    :return: The simplified tables and array mapping input sample IDs to output IDS

    :rtype: tuple

    Note that the samples argument is agnostic with respect to the time of
    the nodes in the input tables. Thus, you may do things like simplify
    to a set of "currently-alive" nodes plus some or all ancient samples by
    including some node IDs from
    :attr:`fwdpy11.DiploidPopulation.ancient_sample_metadata`.

    If the input contains ancient samples, and you wish to include them in the output,
    then you need to include their IDs in the samples argument.

    .. note::

        Due to node ID remapping, the metadata corresponding to nodes becomes a bit more
        difficult to look up.  You need to use the output ID map, the original IDs, and
        the population's metadata containers.

    .. deprecated:: 0.3.0

        Prefer :func:`fwdpp.simplify_tables`


    .. versionchanged:: 0.3.0

        Ancient samples are no longer kept by default

    .. versionchanged:: 0.5.0

        No longer requires a :class:`MutationVector` argument.

    """
    import warnings

    warnings.warn(
        "This function is deprecated and will be removed soon. Please use fwdpy11.simplify_tables instead",
        category=FutureWarning,
    )

    ll_t, idmap = fwdpy11._fwdpy11._simplify(pop, samples)
    return fwdpy11._types.TableCollection(ll_t), idmap


def simplify_tables(
    tables: fwdpy11.TableCollection, samples: Union[List, np.ndarray]
) -> Tuple[fwdpy11.TableCollection, np.ndarray]:
    """
    Simplify a TableCollection.

    :param pop: A table collection.
    :type pop: :class:`fwdpy11.TableCollection`
    :param samples: list of samples
    :type list: list-like or array-like

    :returns: A simplified TableCollection and an array containing remapped sample ids.
    :rtype: tuple

    .. versionadded:: 0.3.0

    """
    ll_t, idmap = fwdpy11._fwdpy11._simplify_tables(tables, samples)

    return fwdpy11._types.TableCollection(ll_t), idmap
