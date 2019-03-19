""" :module PseudomonasDotComScraperTest: Test module for PseudomonasDotComScraper."""

# Import class to be tested.
from GenDBScraper.PseudomonasDotComScraper import PseudomonasDotComScraper
from GenDBScraper.PseudomonasDotComScraper import pdc_query,\
                                                  _dict_to_pdc_query,\
                                                  _simple_get,\
                                                  _pandas_references,\
                                                  _get_bib_from_doi

# Utilities
from TestUtilities.TestUtilities import _remove_test_files

# Alias for generic tests.
TestedClass = PseudomonasDotComScraper

# 3rd party imports
import unittest
import inspect
import os, sys
import pandas
import shutil
import bs4
from bs4 import BeautifulSoup
from pandas import DataFrame
from subprocess import Popen

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

        # Check default attribute values.
        self.assertIsNone(instance._PseudomonasDotComScraper__browser)
        self.assertEqual(instance._PseudomonasDotComScraper__pdc_url, 'https://www.pseudomonas.com')
        self.assertIsInstance(instance._PseudomonasDotComScraper__query, pdc_query)
        self.assertEqual(instance._PseudomonasDotComScraper__query.strain, 'sbw25')

    def test_shaped_constructor_query_namedtuple (self):
        """ Test the shaped constructor (with arguments, query is a namedtuple)."""

        # Instantiate.
        query = pdc_query('sbw25','pflu0916')
        instance = TestedClass(query=query)

        # Check correct class provenance.
        self.assertIsInstance(instance, TestedClass)

        # Check member attribute.
        self.assertEqual(instance._PseudomonasDotComScraper__query, [query])

    def test_shaped_constructor_query_dict (self):
        """ Test the shaped constructor (with arguments, query is a dict)."""

        # Instantiate.
        query = {'strain' : 'sbw25', 'feature' : 'pflu0916'}
        instance = TestedClass(query=query)

        # Check correct class provenance.
        self.assertIsInstance(instance, TestedClass)

        # Check member attribute.
        self.assertEqual(instance._PseudomonasDotComScraper__query, [pdc_query(strain='sbw25', feature='pflu0916')])

    def test_query_setter_named_tuple (self):
        """ Test that the parameter 'query' can be set as a pdc_query."""

        # Parameter to test.
        query = pdc_query(strain='sbw25', feature='pflu0916')

        # Instantiate.
        instance = TestedClass()

        # Set value.
        instance.query = query

        # Check without using property. Query was internally converted to tuple.
        self.assertEqual(instance._PseudomonasDotComScraper__query, [query])

    def test_query_setter_list (self):
        """ Test that the parameter 'query' can be set as a list of pdc_query."""

        # Parameter to test.
        queries = [pdc_query(strain='sbw25', feature='pflu{0:04d}'.format(f)) for f in range(15,18)]

        # Instantiate.
        instance = PseudomonasDotComScraper(query=queries)

        # Check without using property.
        self.assertEqual(instance._PseudomonasDotComScraper__query, queries)

    def test_query_setter_dict (self):
        """ Test that the parameter 'query' can be set as a dict."""

        # Parameter to test.
        query = {'strain':'sbw25', 'feature':'pflu0916'}

        # Instantiate.
        instance = TestedClass()

        # Set value.
        instance.query = query

        # Check without using property.
        self.assertEqual(instance._PseudomonasDotComScraper__query, [pdc_query(strain='sbw25', feature='pflu0916')])

    def test_query_setter_exceptions (self):
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

    def test_query_getter (self):
        """ Test the getter for query. """

        # Parameter to test.
        query = pdc_query(strain='sbw25')

        # Instantiate.
        instance = TestedClass(query)

        # Get value.
        self.assertEqual(instance.query, [query])

    def test_connect (self):
        """ Test the connect and connected logic. """
        # Instantiate the class.
        scraper = PseudomonasDotComScraper(query=pdc_query(strain='sbw25'))
        # connected should be False.
        self.assertFalse(scraper.connected)

        # Connect.
        scraper.connect()

        # Now it should be true.
        self.assertTrue(scraper.connected)

    def test_connect_failure (self):
        """ Test exception upon connection failure. """

        # Instantiate the class.
        scraper = PseudomonasDotComScraper(query=pdc_query(strain='sbw25'))
        # Mess up the url.
        scraper._PseudomonasDotComScraper__pdc_url = "https://some.nonexisting.url"

        # Connect should bail out.
        self.assertRaises(ConnectionError, scraper.connect)

    def test_connected_read_only (self):
        """ Test that the connected status cannot be set. """

        # Instantiate the class.
        scraper = PseudomonasDotComScraper(query=pdc_query(strain='sbw25'))

        with self.assertRaises(AttributeError):
            scraper.connected=True

    def test_query_failure_if_not_connected (self):
        """ Test that the query bails out if not connected to DB."""

        # Instantiate the class.
        scraper = PseudomonasDotComScraper(query=pdc_query(strain='sbw25'))

        # Run the query. This should fail.
        self.assertRaises(RuntimeError, scraper.run_query)

    def test_run_query (self):
        """ Test a method. """

        # Instantiate.
        scraper = PseudomonasDotComScraper(query={'strain':'sbw25', 'feature':'pflu0916'})

        # Connect.
        scraper.connect()

        # Call a method (expected to raise NotImplemented).
        results = scraper.run_query()

        # Check rtype.
        self.assertIsInstance(results, dict)

        # Check query key:
        self.assertIn("sbw25_pflu0916", results.keys())

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

            self.assertIn(heading, results['sbw25_pflu0916'].keys())
            self.assertIsInstance(results['sbw25_pflu0916'][heading], DataFrame)

        # Check content.
        present_indices = results['sbw25_pflu0916']['Gene Feature Overview'].index
        for idx in ['Strain', 'Locus Tag', 'Name', 'Replicon', 'Genomic location']:
            self.assertIn(idx, present_indices)

    def test_dict_to_pdc_query (self):
        """ Test the conversion utility that returns a pdc_query (named_tuple) from a dict. """

        query_dict = {'strain' : 'sbw25', 'feature' : 'pflu0914', 'organism':None}

        query_pdc = _dict_to_pdc_query(**query_dict)

        self.assertIsInstance(query_pdc, pdc_query)
        self.assertEqual(query_pdc.strain, 'sbw25')
        self.assertEqual(query_pdc.feature, 'pflu0914')
        self.assertIsNone(query_pdc.organism)

    def test_results_with_references (self):
        """ Run a query that returns non-empty references. """

        query = pdc_query(strain='UCBPP-PA14', feature='PA14_67210')

        scraper = PseudomonasDotComScraper(query=query)
        scraper.connect()
        results = scraper.run_query()

        print(results['References'])

    def test_pandas_references (self):
        """ Test the references parser function."""

        soup = BeautifulSoup(_simple_get("https://www.pseudomonas.com/feature/show/?id=1661780&view=overview"), 'lxml')

        references = _pandas_references(soup)

        print(references)

    def test_bib_from_doi (self):
        """ Test the get_bib_from_doi utility. """

        doi = '10.1073/pnas.1700286114'

        bib = _get_bib_from_doi(doi)

        expected_keys = ['doi',
                         'first_author',
                         'title',
                         'container',
                         'volume',
                         'page',
                         'date',
                         ]

        present_keys = bib.keys()
        for key in expected_keys:
            self.assertIn(key, present_keys)

    def test_json_io (self):
        """ Test the IO capabilities based on the json format. """

        # Setup the query.
        query = pdc_query(strain='UCBPP-PA14', feature='PA14_67210')

        # Setup the scraper.
        scraper = PseudomonasDotComScraper(query=query)

        # Connect and run the query.
        scraper.connect()
        results = scraper.run_query()

        # Serialize.
        json_path = scraper.to_json(results)
        self._test_files.append(json_path)

        # Check file exists.
        self.assertTrue(os.path.isfile(json_path))

        # Read back in.
        loaded_results = scraper.from_json(json_path)

        # Check entity and keys.
        self.assertIsInstance(loaded_results, dict)
        expected_keys = ["Gene Feature Overview",
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
                        ]
        present_keys = loaded_results.keys()
        for xk in expected_keys:
            self.assertIn(xk, present_keys)
            self.assertIsInstance(loaded_results[xk], pandas.DataFrame)


if __name__ == "__main__":
    unittest.main()
