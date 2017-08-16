# Tests of types in fwdpy11.trait_values

import unittest
import fwdpy11


class testObject_repr(unittest.TestCase):
    def test_SlocusAdditiveTrait(self):
        from fwdpy11.trait_values import SlocusAdditiveTrait
        w = SlocusAdditiveTrait(1.0)
        ww = eval('fwdpy11.' + repr(w))
        self.assertEqual(type(ww), SlocusAdditiveTrait)
        self.assertEqual(ww.scaling, 1.0)

    def test_SlocusMultTrait(self):
        from fwdpy11.trait_values import SlocusMultTrait
        w = SlocusMultTrait(1.0)
        ww = eval('fwdpy11.' + repr(w))
        self.assertEqual(type(ww), SlocusMultTrait)
        self.assertEqual(ww.scaling, 1.0)

    def test_GBR(self):
        from fwdpy11.trait_values import SlocusGBRTrait
        w = SlocusGBRTrait()
        ww = eval('fwdpy11.' + repr(w))
        self.assertEqual(type(ww), SlocusGBRTrait)


if __name__ == "__main__":
    unittest.main()
