import cppimport
cppimport.force_rebuild()
cppimport.set_quiet(False)
snowdrift = cppimport.imp("snowdrift")
test_polymorphism = cppimport.imp("test_polymorphism")
import unittest
import fwdpy11.fitness as fp11w
import re


class testFitnessPolymorphism(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.w = [
            fp11w.SlocusAdditive(),
            fp11w.SlocusMult(),
            snowdrift.SlocusSnowdrift(0.2, -0.2, 1, -2)]
        self.v = test_polymorphism.test_callback_names(self.w)

    def testCallbackNameLen(self):
        self.assertEqual(len(self.v), 3)

    def testReturnValuesSensible(self):
        p1 = re.compile('additive_diploid')
        p2 = re.compile('multiplicative_diploid')
        p3 = re.compile('snowdrift_diploid')
        self.assertEqual(len(p1.findall(self.v[0])), 1)
        self.assertEqual(len(p2.findall(self.v[1])), 1)
        self.assertEqual(len(p3.findall(self.v[2])), 1)

if __name__ == "__main__":
    unittest.main()

