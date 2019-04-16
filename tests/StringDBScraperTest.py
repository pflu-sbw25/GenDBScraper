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
import time

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
        time.sleep(1)

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
        time.sleep(1)

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
        time.sleep(1)

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
        time.sleep(1)

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
        time.sleep(1)

        db = StringDBScraper(query=stringdb_query("216595", ['pflu_4385', 'pflu_0381']))
        db.connect()

        # Get pandas DataFrame of resolved ids.
        resolution = db.resolve_id()

        # Check type and content of returned frame.
        self.assertIsInstance(resolution, pandas.DataFrame)
        self.assertEqual( resolution.index.tolist(), db.query.features)
        self.assertEqual( resolution['preferredName'].loc[db.query.features[0]], db.query.features[0].upper())

    def test_update_features (self):
        """ Test the identifier resolution API with a modified query. """
        time.sleep(1)

        query=stringdb_query("216595", ['pflu_4385', 'pflu_0381'])
        db = StringDBScraper(query=query)
        db.connect()

        # Get pandas DataFrame of resolved ids.
        db.update_features()

        # Check query was updated.
        for i,f in enumerate(query.features):
            self.assertEqual( db.query.features[i], query.features[i].upper())


    def test_resolve_id_query (self):
        """ Test the identifier resolution API with a modified query. """
        time.sleep(1)

        db = StringDBScraper()

        query=stringdb_query("216595", ['pflu_4385', 'pflu_0381'])

        db.connect()

        # Get pandas DataFrame of resolved ids.
        resolution = db.resolve_id(query=query)

        # Check type and content of returned frame.
        self.assertIsInstance(resolution, pandas.DataFrame)
        self.assertEqual( resolution.index.tolist(), db.query.features)
        self.assertEqual( resolution['preferredName'].loc[db.query.features[0]], db.query.features[0].upper())

    @unittest.skip("Don't run me, I'm flodding the server.")
    def test_network_image (self):
        """ Test the network image grabber with various options. """
        time.sleep(1)

        # Setup the API and connect.
        db = StringDBScraper(query=stringdb_query(taxonId="216595", features=['pflu_4385', 'pflu_0381']))
        db.connect()

        # Get network image.
        nodes = [0, 1, 10]
        flavors = ['evidence', 'confidence', 'action']
        formats = ['png', 'image', 'hires_png', 'highres_image', 'svg']


        for wn in nodes:
            for cn in nodes:
                for flavor in flavors:
                    for form in formats:
                        nw_image = db.network_image(show_image=False, image_format=form, white_nodes=wn, color_nodes=cn, flavor=flavor)
                        self.assertTrue(os.path.isfile(nw_image))
                        self._test_files.append(nw_image)
                        time.sleep(1)

    def test_network_image_default (self):
        """ Test the network image grabber. """
        time.sleep(1)

        # Setup the API and connect.
        db = StringDBScraper(query=stringdb_query(taxonId="216595", features=['pflu_4385', 'pflu_0381']))
        db.connect()

        # Get network image.
        nw_image = db.network_image(show_image=False)
        self._test_files.append(nw_image)

        # Check file is written.
        self.assertTrue(os.path.isfile(nw_image))

    def test_network_interactions(self):
        """ Test getting the network interactions as a table. """
        time.sleep(1)

        # Setup the API and connect.
        query = stringdb_query(taxonId="216595", features=['pflu_4385', 'pflu_0381'])
        db = StringDBScraper(query=query)
        db.connect()

        # Get network interactions.
        nw_interactions = db.network_interactions(nodes=10)

        # Check type.
        self.assertIsInstance(nw_interactions, pandas.DataFrame)

        # List of expected column headers.
        column_names = [
                'stringId_A',
                'stringId_B',
                'ncbiTaxonId',
                'score',
                'nscore',
                'fscore',
                'pscore',
                'ascore',
                'escore',
                'dscore',
                'tscore',
                ]

        # Check we have only and all expected columns.
        self.assertEqual(set(column_names), set(nw_interactions.columns))

    def test_interaction_partners(self):
        """ Test getting the interaction partners. """
        time.sleep(1)

        # Setup the API and connect.
        query = stringdb_query(taxonId="216595", features=['pflu_4385', 'pflu_0381'])
        db = StringDBScraper(query=query)
        db.connect()

        # Get network interactions.
        interaction_partners = db.interaction_partners(required_score=300)

        self.assertIsInstance(interaction_partners, pandas.DataFrame)

        column_names = [
                'stringId_A',
                'stringId_B',
                'ncbiTaxonId',
                'score',
                'nscore',
                'fscore',
                'pscore',
                'ascore',
                'escore',
                'dscore',
                'tscore',
                ]

        self.assertEqual(set(column_names), set(interaction_partners.columns))

    def test_similarity_scores(self):
        """ Test getting the similarity scores. """
        time.sleep(1)

        # Setup the API and connect.
        query = stringdb_query(taxonId="216595", features=['pflu_4385', 'pflu_0381'])
        db = StringDBScraper(query=query)
        db.connect()

        # Get network interactions.
        with self.assertRaises(NotImplementedError):
            similarity_scores = db.similarity_scores()
        return

        self.assertIsInstance(similarity_scores, pandas.DataFrame)

        column_names = [
                'stringId_A',
                'stringId_B',
                'bitscore',
                'start_A',
                'end_A',
                'start_B',
                'end_B',
                'size_B',
                ]

        self.assertEqual(set(column_names), set(similarity_scores.columns))

    def test_functional_enrichments(self):
        """ Test getting the functional enrichments. """
        time.sleep(1)

        # Setup the API and connect.
        query = stringdb_query(taxonId="216595", features=['pflu_0438'])
        db = StringDBScraper(query=query)
        db.connect()

        functional_enrichments = db.functional_enrichments()

        self.assertIsInstance(functional_enrichments, pandas.DataFrame)

        column_names = [
                'category',
                'term',
                'number_of_genes',
                'number_of_genes_in_background',
                'ncbiTaxonId',
                'inputGenes',
                'p_value',
                'fdr',
                'description',
                ]

        self.assertEqual(set(column_names), set(functional_enrichments.columns))

    def test_interaction_enrichments(self):
        """ Test getting the functional enrichments. """
        time.sleep(1)

        # Setup the API and connect.
        query = stringdb_query(taxonId="216595", features=['pflu_0438'])
        db = StringDBScraper(query=query)
        db.connect()

        interaction_enrichments = db.interaction_enrichments()

        self.assertIsInstance(interaction_enrichments, pandas.DataFrame)

        column_names = [
                'number_of_nodes',
                'number_of_edges',
                'average_node_degree',
                'local_clustering_coefficient',
                'expected_number_of_edges',
                'p_value',
                ]

        self.assertEqual(set(column_names), set(interaction_enrichments.columns))


if __name__ == "__main__":
    unittest.main()
