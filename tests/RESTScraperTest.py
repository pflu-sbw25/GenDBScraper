""" :module RESTScraperTest: Test module for RESTScraper."""

# Import class to be tested.
from GenDBScraper.RESTScraper import RESTScraper

# Utilities
from TestUtilities.TestUtilities import _remove_test_files
from TestUtilities.TestUtilities import check_keys

# Alias for generic tests.
TestedClass = RESTScraper

import inspect
import os, sys
import shutil
import unittest

class DerivedScraper(RESTScraper):
    """ Derived scraper only needed to run tests on the virtual base class ."""

    def __init__(self):

        base_url = "https://www.rest.api.net"
        # Init base class.
        super(DerivedScraper, self).__init__(base_url)

class RESTScraperTest(unittest.TestCase):
    """ :class: Test class for the RESTScraper """

    @classmethod
    def setUpClass(cls):
        """ Setup the test class. """

        # Setup a list of test files.
        cls._static_test_files = []

    @classmethod
    def tearDownClass(cls):
        """ Tear down the test class. """

        _remove_test_files(cls._static_test_files)

    def setUp (self):
        """ Setup the test instance. """

        # Setup list of test files to be removed immediately after each test method.
        self._test_files = []

    def tearDown (self):
        """ Tear down the test instance. """
        _remove_test_files(self._test_files)

    def test_abc(self):
        """ Test exception upon constructing the abstract base class."""

        # Attempt instantiation.
        self.assertRaises(TypeError, RESTScraper)

    def test_derived_class(self):
        """ Test instantiation of a derived class. """

        dummy = DerivedScraper()
        self.assertIsInstance(dummy, DerivedScraper)
        self.assertIsInstance(dummy, RESTScraper)

        # Check default attribute values.
        self.assertFalse(dummy.connected)
        self.assertEqual(dummy.base_url, "https://www.rest.api.net")

    def test_base_url_exc(self):
        """ Check that base_url is read-only. """

        dummy = DerivedScraper()

        with self.assertRaises(AttributeError):
            dummy.base_url='http://www.google.de'

if __name__ == "__main__":
    unittest.main()

