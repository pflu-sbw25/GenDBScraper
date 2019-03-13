""" :module: Top level Test Suite that contains all Test Classes. """

import unittest
import os, sys

# Import suites to run.
from PseudomonasDotComScraperTest import PseudomonasDotComScraperTest

# Are we running on CI server?
is_travisCI = ("TRAVIS_BUILD_DIR" in list(os.environ.keys())) and (os.environ["TRAVIS_BUILD_DIR"] != "")


# Define the test suite.
def suite():
    suites = [
               unittest.makeSuite(PseudomonasDotComScraperTest, 'test')
             ]

    return unittest.TestSuite(suites)

# Run the top level suite and return a success status code. This enables running an automated git-bisect.
if __name__=="__main__":

    result = unittest.TextTestRunner(verbosity=2).run(suite())

    if result.wasSuccessful():
        print('---> All tests passed. <---')
        sys.exit(0)

    sys.exit(1)

