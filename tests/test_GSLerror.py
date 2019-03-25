import fwdpy11
import gsl_error
import unittest


class testGSLerror(unittest.TestCase):
    def testErrorHandler(self):
        with self.assertRaises(RuntimeError):
            gsl_error.trigger_error()
        try:
            gsl_error.trigger_error()
        except RuntimeError as e:
            self.assertEqual(str(e), "matrix not square")

    def testUnhandledGSLError(self):
        with self.assertRaises(RuntimeError):
            gsl_error.trigger_error_not_handled()


if __name__ == "__main__":
    unittest.main()
