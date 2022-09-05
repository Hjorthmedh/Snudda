import unittest


class MyTestCase(unittest.TestCase):

    def setUp(self):

        # Create config files for the network using ballanddoublestick, and ballanddoublestick_degenerated
        # (Here axons are elongated from 220 to 300 micrometers, and dendrites are shrunk from 200 to 150 micrometers)

        # Create placement files for the WT and degenerated networks, important that they both have the same neurons.

        # Do touch detection both networks

        # Do swap_to_degenerated_morphologies_extended


        pass

    def test_something(self):
        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
