import unittest
import os
from snudda.utils.snudda_path import snudda_parse_path, get_snudda_data, snudda_simplify_path


class TestSnuddaPath(unittest.TestCase):

    def test_parse(self):

        self.assertEqual(snudda_parse_path("$DATA/test", "/abc"), "/abc/test")

        self.assertEqual(snudda_parse_path("$SNUDDA_DATA/test", "abc/def/"),
                         os.path.abspath(os.path.join(os.path.dirname(__file__), "abc/def/test")))

    def test_get(self):

        self.assertEqual(get_snudda_data(),
                         os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "snudda", "data")))

        self.assertEqual(get_snudda_data(snudda_data=os.path.dirname(__file__)),
                         os.path.dirname(__file__))

        # A few more cases to test...

    def test_simplify(self):

        self.assertEqual(snudda_simplify_path("/abc/def", "/abc"), "$SNUDDA_DATA/def")


if __name__ == '__main__':
    unittest.main()
