# Tests of the fwdpy11.fitness module

import unittest
import fwdpy11


class testObject_repr(unittest.TestCase):
    def test_SlocusAdditive(self):
        from fwdpy11.fitness import SlocusAdditive
        w = SlocusAdditive(2.0)
        ww = eval('fwdpy11.' + repr(w))
        self.assertEqual(type(ww), SlocusAdditive)
        self.assertEqual(ww.scaling, 2.0)
        w = SlocusAdditive(1.0)
        ww = eval('fwdpy11.' + repr(w))
        self.assertEqual(type(ww), SlocusAdditive)
        self.assertEqual(ww.scaling, 1.0)

    def test_SlocusMult(self):
        from fwdpy11.fitness import SlocusMult
        w = SlocusMult(1.0)
        ww = eval('fwdpy11.' + repr(w))
        self.assertEqual(type(ww), SlocusMult)


if __name__ == "__main__":
    unittest.main()
