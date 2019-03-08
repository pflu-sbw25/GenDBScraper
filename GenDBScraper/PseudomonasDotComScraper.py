""" :module PseudomonasDotComScraper: Hosting the PseudomonasDotComScraper, an API for the https://www.pseudomonas.com database web interface. """

# Local imports
# import GenDBScraper
# import SeleniumScraper

# 3rd party imports
import selenium

class PseudomonasDotComScraper():
    """ :class: An API for the pseudomonas.com genome database using web scraping technology. """

    # Class constructor
    def __init__(self,
            url=None,
            feature=None,
            ):
        """
        PseudomonasDotComScraper constructor.

        :param url: The URL at which to start scraping. Default: https://www.pseudomonas.com
        :type url: str

        :param feature: The feature to query from the database
        :type feature: str
        """

        # Base class initialization.
        #super(<+ClassName+>).__init__(<+base_class_args+>)

        # Initialize all variables.
        self.__url = None
        self.__feature = None

        # Set attributes via setter.
        self.url = url
        self.feature = feature

    # Attribute accessors
    @property
    def url(self):
        """ Get the initial url.

        :return: The initial URL where scraping begins.
        :rtype:  str

        """

        return self.__url

    @url.setter
    def url(self, val):
        """ Set the initial URL to given value.

        :param val: The value to set.
        :type val: str

        """

        # Checks
        self.__url = val

    @property
    def feature(self):
        """ Get the feature to query.

        :return: The query string.
        :rtype:  str

        """

        return self.__feature

    @feature.setter
    def feature(self, val):
        """ Set the query feature to a given value.

        :param val: The value to set.
        :type val: str

        """

        # Checks
        self.__feature = val



