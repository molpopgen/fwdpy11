import fwdpy11
import cppimport
cppimport.force_rebuild()
gsl = cppimport.imp("gsl_error")
import unittest

class testGSLerror(unittest.TestCase):
    def testErrorHandler(self):
        with self.assertRaises(RuntimeError):
            gsl.trigger_error()

if __name__ == "__main__":
    unittest.main()
