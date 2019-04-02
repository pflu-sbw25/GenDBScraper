""" :module StringDBScraperTest: Test module for StringDBScraper."""

# Import class to be tested.
from GenDBScraper.StringDBScraper import StringDBScraper, stringdb_query
from GenDBScraper.RESTScraper import RESTScraper
from GenDBScraper.Utilities import web_utilities

# Utilities
from TestUtilities.TestUtilities import _remove_test_files
from TestUtilities.TestUtilities import check_keys

# Alias for generic tests.
TestedClass = StringDBScraper

# 3rd party imports
import inspect
import os, sys
import pandas
import re
import shutil
import unittest

class StringDBScraperTest(unittest.TestCase):
    """ :class: Test class for the StringDBScraper """

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

    def test_default_constructor (self):
        """ Test the default class constructor."""

        # Instantiate.
        instance = TestedClass()

        # Check correct class provenance.
        self.assertIsInstance(instance, TestedClass)
        self.assertIsInstance(instance, RESTScraper)

        # Check default attribute values.
        self.assertEqual(instance.base_url, 'http://string-db.org')
        self.assertEqual(instance.query.taxonId, None)
        self.assertEqual(instance.query.features, [])

    def test_shaped_constructor (self):
        """ Test the shaped class constructor."""

        # Instantiate.
        instance = TestedClass(query=stringdb_query(taxonId='216595', features=['pflu5436', 'DnaA']))

        # Check correct class provenance.
        self.assertIsInstance(instance, TestedClass)
        self.assertIsInstance(instance, RESTScraper)

        # Check default attribute values.
        self.assertEqual(instance.base_url, 'http://string-db.org')
        self.assertEqual(instance.query.taxonId, '216595')
        self.assertEqual(instance.query.features, ['pflu5436', 'DnaA'])

    def test_shaped_constructor_dict (self):
        """ Test the shaped class constructor with a query dict."""

        # Instantiate.
        instance = TestedClass(query=dict(taxonId='216595', features=['pflu5436', 'DnaA']))

        # Check correct class provenance.
        self.assertIsInstance(instance, TestedClass)
        self.assertIsInstance(instance, RESTScraper)

        # Check default attribute values.
        self.assertEqual(instance.base_url, 'http://string-db.org')
        self.assertEqual(instance.query.taxonId, '216595')
        self.assertEqual(instance.query.features, ['pflu5436', 'DnaA'])

    def test_connection(self):
        """ Test connecting to the DB. """

        # Instantiate.
        instance = TestedClass(query=stringdb_query("216595", ['pflu4385', 'pflu0325']))

        # We're not connected.
        self.assertFalse(instance.connected)

        # Connect.
        instance.connect()

        # Now we're connected.
        self.assertTrue(instance.connected)

    def test_resolve_id (self):
        """ Test the identifier resolution API. """
        db = StringDBScraper(query=stringdb_query("216595", ['pflu_4385', 'pflu_0381']))
        db.connect()

        # Get pandas DataFrame of resolved ids.
        resolution = db.resolve_id()

        print(resolution)
        # Check type and content of returned frame.
        self.assertIsInstance(resolution, pandas.DataFrame)
        self.assertEqual( resolution.index.tolist(), db.query.features)
        self.assertEqual( resolution['preferredName'].loc[db.query.features[0]], db.query.features[0].upper())

if __name__ == "__main__":
    unittest.main()
