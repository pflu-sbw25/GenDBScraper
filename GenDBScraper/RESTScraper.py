""" :module RESTScraper: Hosting the the abstract base class (abc) for all database "scrapers" using REST-ful APIs."""

from GenDBScraper.Utilities.web_utilities import guarded_get

from abc import ABC, abstractmethod
import logging

# Configure logging.
logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO)

class RESTScraper(ABC):
    """  The abstract base class for all scrapers using RESTful APIs. """

    # The virtual constructor, to be implemented on derived class.
    @abstractmethod
    def __init__(self, base_url):
        """ """
        """ The virtual constructor of the RESTScraper abstract base class. """

        if isinstance(base_url, str):
            self.__base_url = base_url
        else:
            raise TypeError("base_url must be a string.")

        self.__connected = False
        self.__query = None

    @property
    def base_url(self):
        """ Return the base url of the scraper. """
        return self.__base_url

    @base_url.setter
    def base_url(self, value):
        raise AttributeError("base_url is a read-only property.")

    @property
    def connected(self):
        """ Return the connection state. """
        return self.__connected

    @connected.setter
    def connected(self, value):
        raise AttributeError("connected is a read-only property.")

    @property
    @abstractmethod
    def query(self):
        """ Return the query. """
        return self.__query

    @query.setter
    @abstractmethod
    def query(self, value):
        self.__query = value

    def connect(self):
        """ Connect to the database. """

        response = guarded_get(self.base_url)
        self.__connected = True
