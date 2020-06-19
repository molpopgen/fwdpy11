#
# Copyright (C) 2017-2020 Kevin Thornton <krthornt@uci.edu>
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

import unittest

import attr

import fwdpy11.class_decorators


@fwdpy11.class_decorators.attr_class_to_from_dict
@fwdpy11.class_decorators.attr_add_asblack
@attr.s(auto_attribs=True)
class FauxClass(object):
    a: float
    b: float


class TestClassDecorators(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.f = FauxClass(1.0, 3.0)

    def test_black_pretty_print(self):
        try:
            _ = self.f.asblack()
        except:  # NOQA
            self.fail("self.f.asblack() raised unexpected exception")

    def test_to_from_dict(self):
        d = self.f.asdict()
        self.assertEqual(d["a"], self.f.a)
        self.assertEqual(d["b"], self.f.b)
        f = FauxClass.fromdict(d)
        self.assertEqual(f.a, self.f.a)
        self.assertEqual(f.b, self.f.b)


if __name__ == "__main__":
    unittest.main()
