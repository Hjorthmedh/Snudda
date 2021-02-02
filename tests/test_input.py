import unittest
import os


class MyTestCase(unittest.TestCase):

    def setUp(self):

        os.chdir(os.path.dirname(__file__))

        # Setup network so we can test input generation


        pass

    #def test_something(self):
    #    self.assertEqual(True, False)

    # TODO: Test the disabling of input by adding a ! mark in config file before input name

if __name__ == '__main__':
    unittest.main()
