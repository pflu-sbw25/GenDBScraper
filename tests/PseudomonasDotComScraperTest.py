""" :module PseudomonasDotComScraperTest: Test module for PseudomonasDotComScraper."""

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
import numpy
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
    query = pdc_query(strain='sbw25', feature='pflu0916')

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
        self.assertIsInstance(instance._PseudomonasDotComScraper__query[0], pdc_query)
        self.assertEqual(instance._PseudomonasDotComScraper__query[0].strain, 'sbw25')

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
        """ Test running a query. """

        # Instantiate.
        scraper = PseudomonasDotComScraper(query={'strain':'sbw25', 'feature':'pflu0916'})

        # Connect.
        scraper.connect()

        # Call a method (expected to raise NotImplemented).
        scraper.run_query()
        results = scraper.results

        # Check rtype.
        self.assertIsInstance(results, dict)

        # Check query key:
        self.assertIn("sbw25__pflu0916", results.keys())

        # Check keys.
        expected_keys = ["Overview",
                        "Function/Pathways/GO",
                        "Sequences",
                        "Motifs",
                        "Operons",
                        "Transposon Insertions",
                        "Updates",
                        "Orthologs",
                        ]

        check_keys(self, expected_keys, results["sbw25__pflu0916"])

        # Check content.
        present_indices = results['sbw25__pflu0916']['Overview']['Gene Feature Overview'].loc[:,0].values

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

    def test_get_subcellular_localizations (self):
        """ Test the subcellular_localizaton scraping. """

        scraper = setup_scraper_complete()

        soup =  BeautifulSoup(guarded_get("https://www.pseudomonas.com/feature/show/?id=1466562"), 'lxml')
        panels = scraper._get_subcellular_localizations(soup)

        self.assertIsInstance(panels, dict)
        expected_keys = ["Individual Mappings", "Additional evidence"]
        check_keys(self, expected_keys, panels)

    def test_get_overview(self):
        """ Test the scraping of the "Overview" tab on pseudomonas.com."""

        scraper = setup_scraper_complete()

        panels = dict()
        overview_panel = scraper._get_overview("https://www.pseudomonas.com/feature/show/?id=1661780")

        expected_keys = ["Gene Feature Overview",
                         "Cross-References",
                         "Product",
                         "Subcellular Localizations",
                         "Pathogen Association Analysis",
                         "References",
                         ]

        check_keys(self, expected_keys, overview_panel)

        # Check a reference.
        cross_references = overview_panel["Cross-References"]

        refseq = cross_references.loc[cross_references.loc[:,'type']=="RefSeq",'id'].loc[0]

        self.assertEqual(refseq, "YP_793558.1")

    def test_get_sequence(self):
        """ Test the scraping of the "Sequences" tab on pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panel = scraper._get_sequences("https://www.pseudomonas.com/feature/show/?id=1661780" )

        # Check keys.
        self.assertIsInstance(panel, pandas.DataFrame)
        expected_indices = [
                "DNA Sequence Upstream of Gene",
                "DNA Sequence for Gene",
                "DNA Sequence Downstream of Gene",
                "Amino Acid Sequence",
                ]
        indices = list(panel[0])
        for ei in expected_indices:
            self.assertIn(ei, indices)

    def test_get_sequence_fasta(self):
        """ Test the scraping of the "Sequences" tab on pseudomonas.com and seqio parsing"""

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panel = scraper._get_sequences("https://www.pseudomonas.com/feature/show/?id=1661780" )

        # Check gene.
        stream = StringIO(panel.loc[panel[0] == "DNA Sequence for Gene",1].values[0])
        gene_record = SeqIO.read(stream, 'fasta')

        self.assertEqual(gene_record.seq, "ATGCGAATCTCTATTGGTCTATTCATTTTCCTGTTGAGTTTCGGAGTTCCCGCTATGGCCGACAGTAAGCCGTTCATCTGCGTTAATGAGAAAGACCATGCCACCTCTCGATCCCCAGGCCGATGCCTGGTATCGGGAAGCAGTTGCGCTAGCTAAGCCTGACGCCTTGCGTCCTTGGGGACGTATTGTGGACTTATAGTAAGGCAGTTGAGCGTGGGCATTGGAAGGCGATGCATAATTTGGCGAATCTTTATCGCACAGGGTGGCCCGGAGGGGTAGAAAAAGATACGCAGAACATTGGATCTCTATCAAAAGATGATCGATCTGGAGGTGCCCCAAGGGTTCTATGATATGGGAGCAATGATCGGCAATCGTGCAGGGGTCATGAATCCTAACTGACGGGCTTAGTTTTCTTAATAAGGCTGCTAGCCTAGGAAATCCGCCGGCATTAACCGAGCTAGGTAAGCTCTATATATATGTGGCCAAAAAAAGATTTGGGGTTGGCGTATACTCACTGTGCTGCTAGCCAGGGCTATGCGCCGGCTAGTTATGAGTTGGGGGCGTATTACAAGATAGTAGAGCATAATTTCAAAAGCATTGGGTTATTATCAGGCGTCAGTCTCTCAGGGCGGAAAGAGTGCGGCTTTATTTATCTCCGGTGTTTTTGATAAAGCCAGTCCTGATGTCTAGAATGTGGTACGCACCCGATGAGAAATTGCGCAAATTATATGATGGTATTTACGATAAACTTGCCGCTGATCCTGATTTTCGTTTTCCCAACTTGAAAGGACCATCCTCTACCTTCTCACCCGACCCAGGGCTACGATGCAGATCGGCCCGACTGGAAACCGGGGCAGTGA")

    def test_get_functions_pathways_go(self):
        """ Test the scraping of the "Function/Pathways/GO" tab on pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_functions_pathways_go("https://www.pseudomonas.com/feature/show/?id=1661780" )

        # Check keys.
        expected_keys = ["Gene Ontology",
                         "Functional Classifications Manually Assigned by PseudoCAP",
                         "Functional Predictions from Interpro",
                         ]

        check_keys(self, expected_keys, panels)

        # Check the E-values are floats, not strings.
        # Get E-values.
        e_values = panels["Functional Predictions from Interpro"]["E-value"].values
        for ev in e_values:
            self.assertTrue(numpy.issubdtype(ev, numpy.number))

    def test_get_motifs(self):
        """ Test the scraping of the "Motifs" tab on pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_motifs("https://www.pseudomonas.com/feature/show/?id=1661780")

        # Check keys.
        expected_keys = []

        check_keys(self, expected_keys, panels)

    def test_get_operons_missing(self):
        """ Test the scraping of the "operons" tab on pseudomonas.com in the case no operons available."""

        # Get test scraper.
        scraper = setup_scraper_incomplete()

        gene_url = scraper._get_feature_url(scraper.query[0])

        # Get sequences tables.
        panels = scraper._get_operons(gene_url)

        self.assertEqual(panels, dict())

    def test_get_operons(self):
        """ Test the scraping of the "operons" tab on pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_operons("https://www.pseudomonas.com/feature/show/?id=1661780")

        # Check in first operon.
        expected_keys = ['PA14_67230-PA14_67220-PA14_67210-PA14_67200-PA14_67190-PA14_67180']
        check_keys(self, expected_keys, panels)

        # Check in first operon.
        expected_keys = ['Genes', 'Meta', 'References']
        check_keys(self, expected_keys, panels['PA14_67230-PA14_67220-PA14_67210-PA14_67200-PA14_67190-PA14_67180'])

    def test_get_transposon_insertions_missing(self):
        """ Test the scraping of the "Transposon Insertions" tab on pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_incomplete()

        # Get sequences tables.
        panels = scraper._get_transposon_insertions(scraper._get_feature_url(scraper.query[0]))

        # Check keys.
        expected_keys = ['Transposon Insertions in SBW25', 'Transposon Insertions in Orthologs']
        check_keys(self, expected_keys, panels)

    def test_get_transposon_insertions_in_ortholog(self):
        """ Test the scraping of the "Transposon Insertions" tab on pseudomonas.com with TIs in orthologs."""

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_transposon_insertions("https://www.pseudomonas.com/feature/show/?id=1661770")

        # Check keys.
        expected_keys = ['Transposon Insertions in UCBPP-PA14', 'Transposon Insertions in Orthologs']
        check_keys(self, expected_keys, panels)

        ti_orth = panels['Transposon Insertions in Orthologs']
        self.assertEqual(len(ti_orth.columns), 7)
        self.assertEqual(len(ti_orth.index), 2)

        # Mutant IDs should be different
        self.assertEqual(ti_orth.loc[0,'Mutant ID'], 'UWGC: PW9532')
        self.assertEqual(ti_orth.loc[1,'Mutant ID'], 'UWGC: PW9531')

        # Genetic positions should be different
        self.assertEqual(ti_orth.loc[0,'Genomic Position'], '5722847')
        self.assertEqual(ti_orth.loc[1,'Genomic Position'], '5723351')

    def test_get_transposon_insertions(self):
        """ Test the scraping of the "Transposon Insertions" tab on pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_transposon_insertions("https://www.pseudomonas.com/feature/show/?id=1661780")

        # Check keys.
        expected_keys = ['Transposon Insertions in UCBPP-PA14', 'Transposon Insertions in Orthologs']
        check_keys(self, expected_keys, panels)

    def test_get_updates_missing(self):
        """ Test the scraping of the "updates" tab on pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_incomplete()

        # Get sequences tables.
        panels = scraper._get_updates(scraper._get_feature_url(scraper.query[0]))

        # Check keys.
        expected_keys = ["Annotation Updates"]
        check_keys(self, expected_keys, panels)

    def test_get_updates(self):
        """ Test the scraping of the "updates" tab on pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_updates("https://www.pseudomonas.com/feature/show/?id=1661780")

        # Check keys.
        expected_keys = ["Annotation Updates"]
        check_keys(self, expected_keys, panels)

    def test_get_orthologs(self):
        """ Test the scraping of the "Otholog Group Members" tab on pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_orthologs("https://www.pseudomonas.com/feature/show/?id=1661780")

        # Check keys.
        expected_keys = ['Ortholog cluster', 'Ortholog group', 'Ortholog xml']

        check_keys(self, expected_keys, panels)

    def test_get_ortholog_cluster(self):
        """ Test the scraping of the "Otholog Cluster Members" from pseudoluge.pseudomonas.com."""

        # Get test scraper.
        scraper = setup_scraper_incomplete()

        # Get sequences tables.
        panels = scraper._get_orthologs("https://www.pseudomonas.com/feature/show/?id=1459889")

        # Check keys.
        expected_keys = ['Ortholog cluster', 'Ortholog group', 'Ortholog xml']

        check_keys(self, expected_keys, panels)

    def test_get_cross_references(self):
        """ Test the get_cross_references method. """

        # Get test scraper.
        scraper = setup_scraper_complete()

        # Get sequences tables.
        panels = scraper._get_cross_references("https://www.pseudomonas.com/feature/show/?id=1661780")

        urls = panels['url']
        self.assertIn("http://www.ncbi.nlm.nih.gov/protein/YP_793558.1", urls.values)

    def test_results_with_all_tabs (self):
        """ Run a query on a strain with a 'complete' dataset."""

        # Setup the query.
        query = pdc_query(strain='UCBPP-PA14', feature='PA14_67210')

        # Fetch results.
        scraper = PseudomonasDotComScraper(query=query)
        scraper.connect()

        # Setup the query.
        query = pdc_query(strain='UCBPP-PA14', feature='PA14_67210')

        # Fetch results.
        scraper = PseudomonasDotComScraper(query=query)
        scraper.connect()
        scraper.run_query()
        results = scraper.results['UCBPP-PA14__PA14_67210']

        expected_keys = ["Overview",
                         "Sequences",
                         "Function/Pathways/GO",
                         "Motifs",
                         "Operons",
                         "Transposon Insertions",
                         "Updates",
                         "Orthologs",
                         ]

        check_keys(self, expected_keys, results)

        overview = results["Overview"]
        expected_keys = [
                         "Gene Feature Overview",
                         "Cross-References",
                         "Product",
                         "Subcellular Localizations",
                         "Pathogen Association Analysis",
                         "References",
                         ]
        check_keys(self, expected_keys, overview)

        go = results["Function/Pathways/GO"]
        expected_keys = [
                         "Gene Ontology",
                         "Functional Classifications Manually Assigned by PseudoCAP",
                         "Functional Predictions from Interpro",
                        ]
        check_keys(self, expected_keys, go)

    def test_results_with_references (self):
        """ Run a query that returns non-empty references. """

        # Setup the query.
        query = pdc_query(strain='UCBPP-PA14', feature='PA14_67210')

        # Fetch results.
        scraper = PseudomonasDotComScraper(query=query)
        scraper.connect()
        scraper.run_query()
        results = scraper.results['UCBPP-PA14__PA14_67210']["Overview"]['References']

        # Check column names.
        self.assertIn('citation', results.columns)
        self.assertIn('pubmed_url', results.columns)

        # Check first author.
        citation = results.loc[0]['citation']
        rx = re.compile('Allsopp\s')

        self.assertIsNotNone(rx.match(citation))

    def test_pandas_references (self):
        """ Test the references parser function."""

        soup = BeautifulSoup(guarded_get("https://www.pseudomonas.com/feature/show/?id=1661780&view=overview"), 'lxml')

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

    @unittest.expectedFailure
    def test_json_io (self):
        """ Test the IO capabilities based on the json format. """

        # Setup the query.
        query = pdc_query(strain='UCBPP-PA14', feature='PA14_67210')

        # Setup the scraper.
        scraper = PseudomonasDotComScraper(query=query)

        # Connect and run the query.
        scraper.connect()

        scraper.run_query()
        results = scraper.results

        # Serialize.
        json_path = scraper.to_json(results)
        self._test_files.append(json_path)

        # Check file exists.
        self.assertTrue(os.path.isfile(json_path))

        # Read back in.
        loaded_results = scraper.from_json(json_path)

        # Check entity and keys.
        self.assertIsInstance(loaded_results, dict)
        expected_query_key = 'UCBPP-PA14__PA14_67210'
        self.assertIn(expected_query_key, loaded_results.keys())
        self.assertIsInstance(loaded_results[expected_query_key], dict)

        # Check keys in query results.
        query_results = loaded_results['UCBPP-PA14__PA14_67210']
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
        present_keys = query_results.keys()

        for xk in expected_keys:
            self.assertIn(xk, present_keys)
            self.assertIsInstance(query_results[xk], pandas.DataFrame)

if __name__ == "__main__":

    unittest.main()
