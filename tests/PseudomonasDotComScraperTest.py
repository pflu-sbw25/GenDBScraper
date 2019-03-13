""" :module PseudomonasDotComScraperTest: Test module for PseudomonasDotComScraper."""

# Import class to be tested.
from GenDBScraper.PseudomonasDotComScraper import PseudomonasDotComScraper
from GenDBScraper.PseudomonasDotComScraper import pdc_query, _dict_to_pdc_query

# Utilities
from TestUtilities.TestUtilities import _remove_test_files

# Alias for generic tests.
TestedClass = PseudomonasDotComScraper

# 3rd party imports
import unittest
import os
import pandas
import shutil
import bs4
from pandas import DataFrame

class PseudomonasDotComScraperTest(unittest.TestCase):
    """ :class: Test class for the PseudomonasDotComScraper """

    @classmethod
    def setUpClass(cls):
        """ Setup the test class. """

        # Setup a list of test files.
        cls._static_test_files = []

    @classmethod
    def tearDownClass(cls):
        """ Tear down the test class. """

        _remove_test_files(cls._static_test_files)

    def setUp(self):
        """ Setup the test instance. """

        # Setup list of test files to be removed immediately after each test method.
        self._test_files = []

    def tearDown(self):
        """ Tear down the test instance. """
        _remove_test_files(self._test_files)

    def test_default_constructor(self):
        """ Test the default class constructor."""

        # Instantiate.
        instance = TestedClass()

        # Check correct class provenance.
        self.assertIsInstance(instance, TestedClass)

        # Check default attribute values.
        self.assertIsNone(instance._PseudomonasDotComScraper__browser)
        self.assertEqual(instance._PseudomonasDotComScraper__pdc_url, 'https://www.pseudomonas.com')
        self.assertIsInstance(instance._PseudomonasDotComScraper__query, pdc_query)
        self.assertEqual(instance._PseudomonasDotComScraper__query.strain, 'sbw25')

    def test_shaped_constructor_query_namedtuple(self):
        """ Test the shaped constructor (with arguments, query is a namedtuple)."""

        # Instantiate.
        query = pdc_query('sbw25','pflu0916')
        instance = TestedClass(query=query)

        # Check correct class provenance.
        self.assertIsInstance(instance, TestedClass)

        # Check member attribute.
        self.assertEqual(instance._PseudomonasDotComScraper__query, query)

    def test_shaped_constructor_query_dict(self):
        """ Test the shaped constructor (with arguments, query is a dict)."""

        # Instantiate.
        query = {'strain' : 'sbw25', 'feature' : 'pflu0916'}
        instance = TestedClass(query=query)

        # Check correct class provenance.
        self.assertIsInstance(instance, TestedClass)

        # Check member attribute.
        self.assertEqual(instance._PseudomonasDotComScraper__query, pdc_query(strain='sbw25', feature='pflu0916'))

    def test_query_setter_named_tuple(self):
        """ Test that the parameter 'query' can be set as a pdc_query."""

        # Parameter to test.
        query = pdc_query(strain='sbw25', feature='pflu0916')

        # Instantiate.
        instance = TestedClass()

        # Set value.
        instance.query = query

        # Check without using property.
        self.assertEqual(instance._PseudomonasDotComScraper__query, query)

    def test_query_setter_dict(self):
        """ Test that the parameter 'query' can be set as a dict."""

        # Parameter to test.
        query = {'strain':'sbw25', 'feature':'pflu0916'}

        # Instantiate.
        instance = TestedClass()

        # Set value.
        instance.query = query

        # Check without using property.
        self.assertEqual(instance._PseudomonasDotComScraper__query, pdc_query(strain='sbw25', feature='pflu0916'))

    def test_query_setter_exceptions(self):
        """ Test that passing wrong values to query setter raises. """

        # Wrong type
        query = 12

        self.assertRaises(TypeError, TestedClass, query)

        # Wrong keys (giving strain, feature and organism to query.)
        query = pdc_query('sbw25', 'pflu0915', 'pseudomonas')
        self.assertRaises(KeyError, TestedClass, query)

        query = pdc_query('sbw25', organism='pseudomonas')
        self.assertRaises(KeyError, TestedClass, query)

        query = pdc_query(feature='pflu0916', organism='pseudomonas')
        self.assertIsInstance(TestedClass(query), PseudomonasDotComScraper)

        # Wrong typed values
        query = pdc_query(strain=12)
        self.assertRaises(TypeError, TestedClass, query)

    def test_query_getter(self):
        """ Test the getter for query. """

        # Parameter to test.
        query = pdc_query(strain='sbw25')

        # Instantiate.
        instance = TestedClass(query)

        # Get value.
        self.assertEqual(instance.query, query)

    def test_connect(self):
        """ Test the connect and connected logic. """
        # Instantiate the class.
        scraper = PseudomonasDotComScraper(query=pdc_query(strain='sbw25'))
        # connected should be False.
        self.assertFalse(scraper.connected)

        # Connect.
        scraper.connect()

        # Now it should be true.
        self.assertTrue(scraper.connected)

    def test_connect_failure(self):
        """ Test exception upon connection failure. """

        # Instantiate the class.
        scraper = PseudomonasDotComScraper(query=pdc_query(strain='sbw25'))
        # Mess up the url.
        scraper._PseudomonasDotComScraper__pdc_url = "https://some.nonexisting.url"

        # Connect should bail out.
        self.assertRaises(ConnectionError, scraper.connect)

    def test_connected_read_only(self):
        """ Test that the connected status cannot be set. """

        # Instantiate the class.
        scraper = PseudomonasDotComScraper(query=pdc_query(strain='sbw25'))

        with self.assertRaises(AttributeError):
            scraper.connected=True

    def test_query_failure_if_not_connected(self):
        """ Test that the query bails out if not connected to DB."""

        # Instantiate the class.
        scraper = PseudomonasDotComScraper(query=pdc_query(strain='sbw25'))

        # Run the query. This should fail.
        self.assertRaises(RuntimeError, scraper.run_query)

    def test_run_query(self):
        """ Test a method. """

        # Instantiate.
        scraper = PseudomonasDotComScraper(query={'strain':'sbw25', 'feature':'pflu0916'})

        # Connect.
        scraper.connect()

        # Call a method (expected to raise NotImplemented).
        results = scraper.run_query()

        # Check rtype.
        self.assertIsInstance(results, dict)

        # Check keys.
        for heading in ["Gene Feature Overview",
                        "Cross-References",
                        "Product",
                        "Subcellular localization",
                        "Pathogen Association Analysis",
                        "Orthologs/Comparative Genomics",
                        "Interactions",
                        "References",
                        "Gene Ontology",
                        "Functional Classifications Manually Assigned by PseudoCAP",
                        "Functional Predictions from Interpro",
                        ]:
            self.assertIn(heading, results.keys())
            self.assertIsInstance(results[heading], DataFrame)

        # Check content.
        present_indices = results['Gene Feature Overview'].index
        for idx in ['Strain', 'Locus Tag', 'Name', 'Replicon', 'Genomic location']:
            self.assertIn(idx, present_indices)

    def test_dict_to_pdc_query(self):
        """ Test the conversion utility that returns a pdc_query (named_tuple) from a dict. """

        query_dict = {'strain' : 'sbw25', 'feature' : 'pflu0914', 'organism':None}

        query_pdc = _dict_to_pdc_query(**query_dict)

        self.assertIsInstance(query_pdc, pdc_query)
        self.assertEqual(query_pdc.strain, 'sbw25')
        self.assertEqual(query_pdc.feature, 'pflu0914')
        self.assertIsNone(query_pdc.organism)

    def test_results_with_references(self):
        """ Run a query that returns non-empty references. """

        query = pdc_query(strain='UCBPP-PA14', feature='PA14_67210')

        scraper = PseudomonasDotComScraper(query=query)
        scraper.connect()
        results = scraper.run_query()

        print(results['References'])


if __name__ == "__main__":
    unittest.main()
