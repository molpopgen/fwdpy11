from typing import List, Tuple, Union

import fwdpy11._fwdpy11
import fwdpy11._types
import numpy as np


def simplify_tables(
    tables: fwdpy11._types.TableCollection, samples: Union[List, np.ndarray]
) -> Tuple[fwdpy11._types.TableCollection, np.ndarray]:
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
