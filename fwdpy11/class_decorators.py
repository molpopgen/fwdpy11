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

import attr


def fromdict(cls):
    """
    Add a classmethod to create an instance from
    a dictionary.
    """
    def _fromdict(cls, d):
        """Build an instance from a dictionary"""
        return cls(**d)

    cls.fromdict = classmethod(_fromdict)
    return cls


@fromdict
def attr_class_to_from_dict(cls):
    """
    When a class is built using attrs, this decorator
    adds two methods, asdict and fromdict.  The latter
    is a classmethod to create an instance.
    """

    def _asdict(self):
        """Return dict representation"""
        return attr.asdict(self)

    cls.asdict = _asdict
    return cls


@fromdict
def attr_class_to_from_dict_no_recurse(cls):
    """
    When a class is built using attrs, this decorator
    adds two methods, asdict and fromdict.  The latter
    is a classmethod to create an instance.

    This version does not recurse into nested objects.
    """

    def _asdict(self):
        """Return dict representation"""
        return attr.asdict(self, recurse=False)

    cls.asdict = _asdict

    return cls
