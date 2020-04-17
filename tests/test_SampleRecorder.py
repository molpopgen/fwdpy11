import unittest

import numpy as np

import fwdpy11


class testSampleRecorder(unittest.TestCase):
    @classmethod
    def setUp(self):
        self.sr = fwdpy11.SampleRecorder()

    def test_add_sample(self):
        self.sr.add_sample(1)
        self.assertEqual(len(self.sr.samples), 1)
        self.assertEqual(self.sr.samples[0], 1)

    def test_assign(self):
        samples = np.array([1, 2, 3, 4], dtype=np.uint32)
        self.sr.assign(samples)
        self.assertEqual(len(self.sr.samples), 4)
        for i in range(4):
            self.assertEqual(self.sr.samples[i], i + 1)


if __name__ == "__main__":
    unittest.main()
