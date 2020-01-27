"""
:module PseudomonasDotComScraperTest:
    Test module for PseudomonasDotComScraper.
"""

# Import class to be tested.
from GenDBScraper.PseudomonasDotComScraper import PseudomonasDotComScraper
from GenDBScraper.Utilities.web_utilities import guarded_get
from GenDBScraper.PseudomonasDotComScraper import pdc_query,\
                                                  _dict_to_pdc_query,\
                                                  _pandas_references,\
                                                  _get_bib_from_doi

# Utilities
from TestUtilities.TestUtilities import _remove_test_files
from TestUtilities.TestUtilities import check_keys
# 3rd party imports
from bs4 import BeautifulSoup
import os
import pandas
import re
import unittest
from io import StringIO
from Bio import SeqIO

# Alias for generic tests.
TestedClass = PseudomonasDotComScraper


def setup_scraper_complete():
    """ Construct a default scraper for testing. """

    # Setup the query.
    query = pdc_query(strain='UCBPP-PA14', feature='PA14_67210')

    # Fetch results.
    scraper = PseudomonasDotComScraper(query=query)
    scraper.connect()

    return scraper


def setup_scraper_incomplete():
    """ Construct a default scraper for testing. """

    # Setup the query.
    query = pdc_query(strain='SBW25', feature='PFLU0916')

    # Fetch results.
    scraper = PseudomonasDotComScraper(query=query)
    scraper.connect()

    return scraper


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

        # Setup list of test files to be removed immediately
        # after each test method.
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
        self.assertEqual(instance._PseudomonasDotComScraper__pdc_url,
                         'https://www.pseudomonas.com'
                         )
        self.assertIsInstance(instance._PseudomonasDotComScraper__query,
                              pdc_query
                              )
        self.assertEqual(instance._PseudomonasDotComScraper__query.strain,
                         'sbw25'
                         )

    def test_shaped_constructor_query_namedtuple(self):
        """
        Test the shaped constructor (with arguments, query is a namedtuple).
        """

        # Instantiate.
        query = pdc_query('SBW25', 'PFLU0916')
        instance = TestedClass(query=query)

        # Check correct class provenance.
        self.assertIsInstance(instance, TestedClass)

        # Check member attribute.
        self.assertEqual(instance._PseudomonasDotComScraper__query, query)

    def test_shaped_constructor_query_dict(self):
        """ Test the shaped constructor (with arguments, query is a dict)."""

        # Instantiate.
        query = {'strain': 'SBW25', 'feature': 'PFLU0916'}
        instance = TestedClass(query=query)

        # Check correct class provenance.
        self.assertIsInstance(instance, TestedClass)

        # Check member attribute.
        self.assertEqual(instance._PseudomonasDotComScraper__query,
                             pdc_query(strain='SBW25', feature='PFLU0916'),
                        )

    def test_query_setter_namedtuple(self):
        """ Test that the parameter 'query' can be set as a pdc_query."""

        # Parameter to test.
        query = pdc_query(strain='SBW25', feature='PFLU0916')

        # Instantiate.
        instance = TestedClass()

        # Set value.
        instance.query = query

        # Check without using property.
        # Query was internally converted to tuple.
        self.assertEqual(instance._PseudomonasDotComScraper__query, query)

    def test_query_setter_list(self):
        """
        Test that the parameter 'query' can be set as a list of pdc_query.
        """

        # Parameter to test.
        queries = [pdc_query(strain='SBW25',
                             feature='PFLU{0:04d}'.format(f))
                   for f in range(15, 18)
                   ]

        # Instantiate.
        self.assertRaises(TypeError, PseudomonasDotComScraper, queries)

    def test_query_setter_dict(self):
        """ Test that the parameter 'query' can be set as a dict."""

        # Parameter to test.
        query = {'strain': 'SBW25', 'feature': 'PFLU0916'}

        # Instantiate.
        instance = TestedClass()

        # Set value.
        instance.query = query

        # Check without using property.
        self.assertEqual(instance._PseudomonasDotComScraper__query,
                             pdc_query(strain='SBW25', feature='PFLU0916')
                         )

    def test_query_setter_exceptions(self):
        """ Test that passing wrong values to query setter raises. """

        # Wrong type
        query = 12

        self.assertRaises(TypeError, TestedClass, query)

        # Wrong keys (giving strain, feature and organism to query.)
        query = pdc_query('SBW25', 'PFLU0915', 'pseudomonas')
        self.assertRaises(KeyError, TestedClass, query)

        query = pdc_query('SBW25', organism='pseudomonas')
        self.assertRaises(KeyError, TestedClass, query)

        query = pdc_query(feature='PFLU0916', organism='pseudomonas')
        self.assertIsInstance(TestedClass(query), PseudomonasDotComScraper)

        # Wrong typed values
        query = pdc_query(strain=12)
        self.assertRaises(TypeError, TestedClass, query)

    def test_query_getter(self):
        """ Test the getter for query. """

        # Parameter to test.
        query = pdc_query(strain='SBW25')

        # Instantiate.
        instance = TestedClass(query)

        # Get value.
        self.assertEqual(instance.query, query)

    def test_connect(self):
        """ Test the connect and connected logic. """
        # Instantiate the class.
        scraper = PseudomonasDotComScraper(query=pdc_query(strain='SBW25'))
        # connected should be False.
        self.assertFalse(scraper.connected)

        # Connect.
        scraper.connect()

        # Now it should be true.
        self.assertTrue(scraper.connected)

    def test_connect_failure(self):
        """ Test exception upon connection failure. """

        # Instantiate the class.
        scraper = PseudomonasDotComScraper(query=pdc_query(strain='SBW25'))
        # Mess up the url.
        scraper._PseudomonasDotComScraper__pdc_url = \
            "https://some.nonexisting.url"

        # Connect should bail out.
        self.assertRaises(ConnectionError, scraper.connect)

    def test_connected_read_only(self):
        """ Test that the connected status cannot be set. """

        # Instantiate the class.
        scraper = PseudomonasDotComScraper(query=pdc_query(strain='SBW25'))

        with self.assertRaises(AttributeError):
            scraper.connected = True

    def test_query_failure_if_not_connected(self):
        """ Test that the query bails out if not connected to DB."""

        # Instantiate the class.
        scraper = PseudomonasDotComScraper(query=pdc_query(strain='SBW25'))

        # Run the query. This should fail.
        self.assertRaises(RuntimeError, scraper.run_query)

    def test_run_query(self):
        """ Test running a query. """

        # Instantiate.
        scraper = PseudomonasDotComScraper(
            query={'strain': 'SBW25',
                   'feature': 'PFLU0916',
                   }
        )

        # Connect.
        scraper.connect()

        # Call a method (expected to raise NotImplemented).
        scraper.run_query()
        results = scraper.results

        # Check rtype.
        self.assertIsInstance(results, dict)

        # Check query key:
        self.assertIn("SBW25__PFLU0916", results.keys())

        # Check keys.
        expected_keys = ["Overview",
                         "Function/Pathways/GO",
                         "Operons",
                         "Transposon Insertions",
                         ]

        check_keys(self, expected_keys, results["SBW25__PFLU0916"])

    def test_dict_to_pdc_query(self):
        """
        Test the conversion utility that returns a pdc_query (named_tuple)
        from a dict.
        """

        query_dict = {'strain': 'SBW25',
                      'feature': 'PFLU0914',
                      'organism': None}

        query_pdc = _dict_to_pdc_query(**query_dict)

        self.assertIsInstance(query_pdc, pdc_query)
        self.assertEqual(query_pdc.strain, 'SBW25')
        self.assertEqual(query_pdc.feature, 'PFLU0914')
        self.assertIsNone(query_pdc.organism)

    def test_get_subcellular_localizations(self):
        """ Test the subcellular_localizaton scraping. """

        scraper = setup_scraper_complete()

        soup = BeautifulSoup(guarded_get(
            "https://www.pseudomonas.com/feature/show/?id=1466562"),
            'lxml'
        )
        panels = scraper._get_subcellular_localizations(soup)

        self.assertIsInstance(panels, dict)
        expected_keys = ["Individual Mappings", "Additional evidence"]
        check_keys(self, expected_keys, panels)

    def test_get_overview(self):
        """ Test the scraping of the "Overview" tab on pseudomonas.com."""

        scraper = setup_scraper_complete()

        overview_panel = scraper._get_overview(
            "https://www.pseudomonas.com/feature/show/?id=1661780"
        )

        expected_keys = ['Strain',
                         'Locus_Tag',
                         'Name',
                         'Replicon',
                         'Transposon_Mutants',
                         'Start',
                         'Stop',
                         'Strand',
                         'RefSeq',
                         'GI',
                         'Entrez',
                         'NCBILocusTag',
                         'Feature_Type',
                         'Coding_Frame',
                         'ProductName',
                         'Synonyms',
                         'Evidence_for_Translation',
                         'Charge_(pH_7)',
                         'Kyte-Doolittle_Hydrophobicity_Value',
                         'Molecular_Weight_(kDa)',
                         'Isoelectric_Point_(pI)',
                         'pubmed_references'
                         ]

        check_keys(self, expected_keys, overview_panel)


        refseq = overview_panel['RefSeq']
        self.assertEqual(refseq, "YP_793558.1")

        product = overview_panel["ProductName"]
        self.assertEqual(product,
                         """type VI secretion lipase immunity protein, Tli5b1 Product Name Confidence: Class 1"""
                         )

    def test_get_sequence(self):
        """ Test the scraping of the "Sequences" tab on pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panel = scraper._get_sequences(
            "https://www.pseudomonas.com/feature/show/?id=1661780"
        )

        # Check keys.
        self.assertIsInstance(panel, dict)
        expected_indices = [
                "DNA Sequence Upstream of Gene",
                "DNA Sequence for Gene",
                "DNA Sequence Downstream of Gene",
                "Amino Acid Sequence",
                ]
        check_keys(self, expected_indices, panel)

    def test_get_sequence_fasta(self):
        """
        Test the scraping of the "Sequences" tab on pseudomonas.com
        and seqio parsing
        """

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panel = scraper._get_sequences(
            "https://www.pseudomonas.com/feature/show/?id=1661780"
        )

        # Check gene.
        stream = StringIO(panel["DNA Sequence for Gene"])
        gene_record = SeqIO.read(stream, 'fasta')

        self.assertEqual(gene_record.seq,
                         "ATGCGAATCTCTATTGGTCTATTCATTTTCCTGTTGAGTTTCGGAGTTCCCGCTATGGCCGACAGTAAGCCGTTCATCTGCGTTAATGAGAAAGACCATGCCACCTCTCGATCCCCAGGCCGATGCCTGGTATCGGGAAGCAGTTGCGCTAGCTAAGCCTGACGCCTTGCGTCCTTGGGGACGTATTGTGGACTTATAGTAAGGCAGTTGAGCGTGGGCATTGGAAGGCGATGCATAATTTGGCGAATCTTTATCGCACAGGGTGGCCCGGAGGGGTAGAAAAAGATACGCAGAACATTGGATCTCTATCAAAAGATGATCGATCTGGAGGTGCCCCAAGGGTTCTATGATATGGGAGCAATGATCGGCAATCGTGCAGGGGTCATGAATCCTAACTGACGGGCTTAGTTTTCTTAATAAGGCTGCTAGCCTAGGAAATCCGCCGGCATTAACCGAGCTAGGTAAGCTCTATATATATGTGGCCAAAAAAAGATTTGGGGTTGGCGTATACTCACTGTGCTGCTAGCCAGGGCTATGCGCCGGCTAGTTATGAGTTGGGGGCGTATTACAAGATAGTAGAGCATAATTTCAAAAGCATTGGGTTATTATCAGGCGTCAGTCTCTCAGGGCGGAAAGAGTGCGGCTTTATTTATCTCCGGTGTTTTTGATAAAGCCAGTCCTGATGTCTAGAATGTGGTACGCACCCGATGAGAAATTGCGCAAATTATATGATGGTATTTACGATAAACTTGCCGCTGATCCTGATTTTCGTTTTCCCAACTTGAAAGGACCATCCTCTACCTTCTCACCCGACCCAGGGCTACGATGCAGATCGGCCCGACTGGAAACCGGGGCAGTGA",
                         )

    def test_get_functions_pathways_go(self):
        """
        Test the scraping of the "Function/Pathways/GO" tab on pseudomonas.com.
        """

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_functions_pathways_go(
            "https://www.pseudomonas.com/feature/show/?id=1661780"
        )

        # Check keys.
        expected_keys = ["Gene Ontology",
                         "Functional Predictions from Interpro",
                         ]

        check_keys(self, expected_keys, panels)

    def test_get_motifs(self):
        """ Test the scraping of the "Motifs" tab on pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_motifs(
            "https://www.pseudomonas.com/feature/show/?id=1661780"
        )

        # Check keys.
        expected_keys = []

        check_keys(self, expected_keys, panels)

    def test_get_operons_missing(self):
        """
        Test the scraping of the "operons" tab on pseudomonas.com in the case
        no operons available.
        """

        # Get test scraper.
        scraper = setup_scraper_incomplete()

        gene_url = scraper._get_feature_url(scraper.query)

        # Get sequences tables.
        panels = scraper._get_operons(gene_url)

        self.assertEqual(panels, dict())

    def test_get_operons(self):
        """ Test the scraping of the "operons" tab on pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_operons(
            "https://www.pseudomonas.com/feature/show/?id=1661780"
        )

        # Check in first operon.
        expected_keys = [
            'PA14_67230-PA14_67220-PA14_67210-PA14_67200-PA14_67190-PA14_67180'
        ]
        check_keys(self, expected_keys, panels)

        # Check in first operon.
        expected_keys = ['Genes', 'Meta', 'References']
        check_keys(self,
                   expected_keys,
                   panels[
                       'PA14_67230-PA14_67220-PA14_67210-PA14_67200-PA14_67190-PA14_67180'
                   ]
                   )

    def test_get_transposon_insertions_missing(self):
        """
        Test scraping of the "Transposon Insertions" tab on pseudomonas.com.
        """

        # Get test scraper.
        scraper = setup_scraper_incomplete()

        # Get sequences tables.
        panels = scraper._get_transposon_insertions(
            scraper._get_feature_url(scraper.query
                                     )
        )

        self.assertEqual(panels, [])

    def test_get_transposon_insertions(self):
        """
        Test the scraping of the "Transposon Insertions" tab
        on pseudomonas.com with TIs in orthologs.
        """

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_transposon_insertions(
            "https://www.pseudomonas.com/feature/show/?id=1661770"
        )

        self.assertEqual(len(panels), 3)

        # Mutant IDs should be different
        self.assertEqual(panels[0]['Mutant ID'], 'UWGC: PW9532')
        self.assertEqual(panels[1]['Mutant ID'], 'UWGC: PW9531')

        # Genetic positions should be different
        self.assertEqual(panels[0]['Genomic Position'], '5722847')
        self.assertEqual(panels[1]['Genomic Position'], '5723351')

    def test_get_updates_missing(self):
        """ Test the scraping of the "updates" tab on pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_incomplete()

        # Get sequences tables.
        panels = scraper._get_updates(
            scraper._get_feature_url(scraper.query)
        )

        # Check keys.
        expected_keys = ["Annotation Updates"]
        check_keys(self, expected_keys, panels)

    def test_get_updates(self):
        """ Test the scraping of the "updates" tab on pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_updates(
            "https://www.pseudomonas.com/feature/show/?id=1661780"
        )

        # Check keys.
        expected_keys = ["Annotation Updates"]
        check_keys(self, expected_keys, panels)

    def test_get_orthologs(self):
        """
        Test the scraping of the "Ortholog Group Members" tab
        on pseudomonas.com.
        """

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_orthologs(
            "https://www.pseudomonas.com/feature/show/?id=1661780"
        )

        # Check keys.
        expected_keys = ['Ortholog cluster', 'Ortholog group', 'Ortholog xml']

        check_keys(self, expected_keys, panels)

    def test_get_ortholog_cluster(self):
        """ Test the scraping of the "Otholog Cluster Members"
        from pseudoluge.pseudomonas.com.
        """

        # Get test scraper.
        scraper = setup_scraper_incomplete()

        # Get sequences tables.
        panels = scraper._get_orthologs(
            "https://www.pseudomonas.com/feature/show/?id=1459889"
        )

        # Check keys.
        expected_keys = ['Ortholog cluster', 'Ortholog group', 'Ortholog xml']

        check_keys(self, expected_keys, panels)

    def test_get_cross_references(self):
        """ Test the get_cross_references method. """

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_cross_references(
            "https://www.pseudomonas.com/feature/show/?id=1661780"
        )


if __name__ == "__main__":

    unittest.main()
