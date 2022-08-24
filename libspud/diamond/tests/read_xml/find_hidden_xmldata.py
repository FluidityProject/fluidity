#!/usr/bin/env python3
import unittest

import diamond.schema as schema


class FindHiddenModule(unittest.TestCase):
    def testNoMagicComments(self):
        schema.find_hidden_xmldata(self.assertEqual(self.schema.tree.getroot(), None))


if __name__ == "__main__":
    unittest.main()
