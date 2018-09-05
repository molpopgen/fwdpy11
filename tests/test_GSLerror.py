import fwdpy11
import cppimport
cppimport.force_rebuild()
gsl = cppimport.imp("gsl_error")
import unittest


class testGSLerror(unittest.TestCase):
    def testErrorHandler(self):
        with self.assertRaises(RuntimeError):
            gsl.trigger_error()
        try:
            gsl.trigger_error()
        except RuntimeError as e:
            self.assertEqual(str(e), "matrix not square")

    def testUnhandledGSLError(self):
        with self.assertRaises(RuntimeError):
            gsl.trigger_error_not_handled()


if __name__ == "__main__":
    unittest.main()
