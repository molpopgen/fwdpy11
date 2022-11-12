import unittest

import fwdpy11
import gsl_error


class testGSLerror(unittest.TestCase):
    def testErrorHandler(self):
        with self.assertRaises(RuntimeError):
            gsl_error.trigger_error_handle_locally()
        try:
            gsl_error.trigger_error_handle_locally()
        except RuntimeError as e:
            self.assertEqual(str(e), "matrix not square")

    def testUnhandledGSLError(self):
        with self.assertRaises(fwdpy11.GSLError):
            gsl_error.trigger_error()


if __name__ == "__main__":
    unittest.main()
