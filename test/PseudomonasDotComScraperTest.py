""" :module PseudomonasDotComScraperTest: Test module for PseudomonasDotComScraper."""

# Import class to be tested.
from GenDBScraper.PseudomonasDotComScraper import PseudomonasDotComScraper as TestedClass

# 3rd party imports
import unittest
import os
import shutil

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

        # Check further ancestry.

        # Check default attribute values.
        self.assertEqual(instance.url, 'https://www.pseudomonas.com')
        self.assertEqual(instance.feature, '')

    def test_shaped_constructor(self):
        """ Test the shaped constructor (with arguments)."""

        # Instantiate.
        instance = TestedClass(url='https://www.pseudomonas.com/index.html')

        # Check correct class provenance.
        self.assertIsInstance(instance, TestedClass)

    def test_url_setter(self):
        """ Test that the parameter 'url' can be set."""

        # Parameter to test.
        url = 'https://www.pseudomonas.com/index.html'

        # Instantiate.
        instance = TestedClass()

        # Set value.
        instance.url = url

        # Check without using property.
        self.assertEqual(instance._PseudomonasDotComScraper__url, url)

    def test_url_setter_exceptions(self):
        """ Test that passing wrongly typed values to url setter raises. """

        url = 12
        instance = TestedClass()

        self.assertRaises(TypeError, instance.url)

    def test_url_getter(self):
        """ Test the getter for url. """

        # Parameter to test.
        url = 'https://www.nonsense.url'

        # Instantiate.
        instance = TestedClass(url)

        # Get value.
        self.assertEqual(instance.url, url)

    def test_feature_setter(self):
        """ Test that the parameter 'feature' can be set."""

        # Parameter to test.
        feature = 'pflu0084'

        # Instantiate.
        instance = TestedClass()

        # Set value.
        instance.feature = feature

        # Check without using property.
        self.assertEqual(instance._PseudomonasDotComScraper__feature, feature)

    def test_feature_setter_exceptions(self):
        """ Test that passing wrongly typed values to feature setter raises. """

        feature = ['pflu0084', 'pflu0916']
        instance = TestedClass()

        self.assertRaises(TypeError, instance.feature)

    def test_feature_getter(self):
        """ Test the getter for feature. """

        # Parameter to test.
        feature = 'pflu0916'

        # Instantiate.
        instance = TestedClass(feature=feature)

        # Get value.
        self.assertEqual(instance.feature, feature)

    @unittest.expectedFailure
    def test_method(self):
        """ Test a method. """

        # Instantiate.
        instance = TestClass()

        # Call a method (expected to raise NotImplemented).
        instance.method()

def _remove_test_files(files):
    """ """
    """ Remove all files and directories listed.

    :param files: Files and directories to remove.
    :type files: list

    """

    # Loop over files
    for f in files:
        # Check if file.
        if os.is_file(f):
            os.remove(f)
        # Check if dir.
        elif os.is_dir(f):
            shutil.rmtree(f)

if __name__ == "__main__":
    unittest.main()
