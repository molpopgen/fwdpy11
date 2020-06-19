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


def _fromdict(cls):
    """
    Add a classmethod to create an instance from
    a dictionary.
    """

    def _fromdict(cls, d):
        """Build an instance from a dictionary"""
        return cls(**d)

    cls.fromdict = classmethod(_fromdict)
    return cls


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
    return _fromdict(cls)


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

    return _fromdict(cls)


def _add_getstate(cls):
    """
    Add __getstate__ via the __dict__
    """

    def getstate(self):
        return self.asdict()

    cls.__getstate__ = getstate

    return cls


def attr_class_pickle(cls):
    """
    Simplistic pickling based
    on the __dict__
    """

    def setstate(self, d):
        self.__dict__.update(d)

    cls.__setstate__ = setstate

    return _add_getstate(cls)


def attr_class_pickle_with_super(cls):
    """
    Add simplistic pickling to a class
    that inherits from another.

    This is tied to our coding stanard,
    where a Python class inheriting from
    a C++ class must use common kwargs
    for both.

    In some cases, the C++ class requires
    more nuanced initialization, so we
    cannot use this decorator.
    """

    def setstate(self, d):
        self.__dict__.update(d)
        super(cls, self).__init__(**d)

    cls.__setstate__ = setstate
    return _add_getstate(cls)


def attr_add_asblack(cls):
    """
    The default __repr__ from attrs isn't readable
    for complex classes.  This adds a method
    for pretty-printing out the class using black's
    formatting rules
    """

    def _asblack(self):
        """
        Return a string representation formatted with black
        """
        import black

        # The try/except is to handle black's changing API
        try:
            return black.format_str(str(self), mode=black.Mode())
        except AttributeError:
            return black.format_str(str(self), mode=black.FileMode())

    cls.asblack = _asblack

    return cls
