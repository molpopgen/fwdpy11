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

from fwdpy11.class_decorators import (attr_class_to_from_dict,
                                      attr_class_to_from_dict_no_recurse)

_common_attr_attribs = {
    "kw_only": True,
    "frozen": True,
    "auto_attribs": True,
    "repr_ns": "fwdpy11.demographic_models",
}


@attr_class_to_from_dict_no_recurse
@attr.s(**_common_attr_attribs)
class DemographicModelDetails(object):
    model: object
    name: str
    source: typing.Dict
    parameters: object
    citation: typing.Optional[typing.Dict]
    metadata: typing.Optional[object] = None

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(d)


@attr_class_to_from_dict
@attr.s(**_common_attr_attribs)
class DemographicModelCitation(object):
    DOI: object
    full_citation: object
    metadata: object

    def __getstate__(self):
        return self.asdict()

    def __setstate__(self, d):
        self.__dict__.update(d)
