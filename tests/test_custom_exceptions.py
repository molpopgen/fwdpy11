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

import unittest

import fwdpy11


class TestFwdppExceptions(unittest.TestCase):
    def test_tables_error(self):
        with self.assertRaises(fwdpy11.TablesError):
            raise fwdpy11.TablesError("this is a tables error")

    def test_tables_error_from_cpp(self):
        import fwdpy11._fwdpy11
        with self.assertRaises(fwdpy11.TablesError):
            fwdpy11._fwdpy11._throw_TablesError()

    def test_samples_error(self):
        with self.assertRaises(fwdpy11.SamplesError):
            raise fwdpy11.SamplesError("this is a tables error")

    def test_samples_error_from_cpp(self):
        import fwdpy11._fwdpy11
        with self.assertRaises(fwdpy11.SamplesError):
            fwdpy11._fwdpy11._throw_SamplesError()


if __name__ == "__main__":
    unittest.main()
