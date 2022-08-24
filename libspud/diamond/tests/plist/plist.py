#!/usr/bin/env python3
import unittest

import diamond.plist as plist


class PyListModule(unittest.TestCase):

    """#1: A simple list of integers, with cardinality ''. (One element only)."""

    def testSimple(self):
        _list = plist.List(int)
        self.assertEqual(_list.__str__(), "list of <type 'int'> of cardinality: ")
        self.assertEqual(_list.__repr__(), "list of <type 'int'> of cardinality: ")
        self.assertEqual(_list("0"), "0")

    """ #2: A simple list of integers, with cardinality '+'. """

    def testOneOrMore(self):
        _list = plist.List(int, "+")
        self.assertEqual(_list.__str__(), "list of <type 'int'> of cardinality: +")
        self.assertEqual(_list.__repr__(), "list of <type 'int'> of cardinality: +")
        self.assertEqual(_list("3,4,5"), "3 4 5")

    """ #3: A list of two strings, with cardinality 2. """

    def testTwoStrings(self):
        _list = plist.List(str, "2")
        self.assertEqual(_list.__str__(), "list of <type 'str'> of cardinality: 2")
        self.assertEqual(_list.__repr__(), "list of <type 'str'> of cardinality: 2")
        self.assertEqual(_list("first second"), "first second")

    """ #4: A list of none type, which should throw an non-callable exception when \
called. """

    def testNoneType(self):
        _list = plist.List(None)
        try:
            _list("3,4,5")
            self.fail()
        except Exception:
            pass


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(PyListModule)
    unittest.TextTestRunner(verbosity=3).run(suite)
