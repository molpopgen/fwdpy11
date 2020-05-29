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

import typing

import attr

from fwdpy11.class_decorators import (attr_add_asblack, attr_class_pickle,
                                      attr_class_to_from_dict,
                                      attr_class_to_from_dict_no_recurse)

_common_attr_attribs = {
    "kw_only": True,
    "frozen": True,
    "auto_attribs": True,
    "repr_ns": "fwdpy11.demographic_models",
}


@attr_add_asblack
@attr_class_pickle
@attr_class_to_from_dict_no_recurse
@attr.s(**_common_attr_attribs)
class DemographicModelDetails(object):
    """
    Stores rich information about a demographic model.
    Instances of this class get returned by functions
    generating pre-calculated models.

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param model: The demographic model parameters
    :type model: object
    :param name: The name of the model
    :type name: str
    :param source: The source of the model
    :type source: dict
    :param parameters: The parameters used to generate ``model``
    :type parameters: object
    :param citation: Citation information
    :type citation: dict or fwdpy11.demographic_models.DemographicModelCitation or None
    :param metadata: Optional field for additional info
    :type metadata: object

    .. versionadded:: 0.8.0
    """

    model: object
    name: str
    source: typing.Dict
    parameters: object
    citation: typing.Optional[typing.Dict]
    metadata: typing.Optional[object] = None


@attr_add_asblack
@attr_class_pickle
@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class DemographicModelCitation(object):
    """
    Citation information for a demographic model

    This class has the following attributes, whose names
    are also ``kwargs`` for intitialization.  The attribute names
    also determine the order of positional arguments:

    :param DOI: The Digital Object Identifier
    :param full_citation: Something string-like giving the full citation.
    :param metadata: Any additional information that may be needed.


    .. versionadded:: 0.8.0
    """

    DOI: object
    full_citation: object
    metadata: object
