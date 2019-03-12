""" :module PseudomonasDotComScraper: Hosting the PseudomonasDotComScraper, an API for the https://www.pseudomonas.com database web interface. """

# Local imports
# import GenDBScraper
# import SeleniumScraper

# 3rd party imports
from requests import get
from requests.exceptions import RequestException
from contextlib import closing
from bs4 import BeautifulSoup
import bs4
import re

from collections import namedtuple

pdc_query = namedtuple('pdc_query',
        field_names=('strain', 'feature', 'organism'),
        defaults=(None, None, None),
        )

class PseudomonasDotComScraper():
    """ :class: An API for the pseudomonas.com genome database using web scraping technology. """

    # Class constructor
    def __init__(self,
            query=None,
            feature=None,
            ):
        """
        PseudomonasDotComScraper constructor.

        :param query: The query to submit to the database.
        :type query: (pdc_query || dict)

        :example: scraper = PseudomonasDotComScraper(query={'strain' : 'sbw25', 'feature' : 'pflu0916'})
        :example: scraper = PseudomonasDotComScraper(query=pdc_query(strain='sbw25', feature='pflu0916'))

        """

        # Base class initialization.
        #super(<+ClassName+>).__init__(<+base_class_args+>)

        # Initialize all variables.
        self.__query = None
        self.__pdc_url = 'https://www.pseudomonas.com'
        self.__browser = None

        # Set attributes via setter.
        self.query = query

    # Attribute accessors
    @property
    def query(self):
        """ Get the query.

        :return: The strings to query the database.
        :rtype:  str

        """

        return self.__query

    @query.setter
    def query(self, val):
        """ Set query.

        :param val: The value to set.
        :type val: str

        """

        # Checks
        if val is None:
            val = pdc_query(strain='sbw25')

        if not isinstance(val, (dict, pdc_query)):
            raise TypeError("The parameter 'query' must be a dict or pdc_query. Examples: query={'strain' : 'sbw25', 'feature'='pflu0916'}; query=pdc_query(strain='sbw25', feature='pflu0916').")

        # Check keys if dict.
        if isinstance(val, dict):
            # Only these are acceptable query keywords.
            accepted_keys = ('strain', 'feature', 'organism')
            present_keys = val.keys()
            for k in present_keys:
                if not k in accepted_keys:
                    raise KeyError("Only 'strain', 'feature', and 'organism' are acceptable keys.)")

            # Complete keywords.
            if not 'strain' in val.keys():
                val['strain'] = None
            if not 'feature' in val.keys():
                val['feature'] = None
            if not 'organism' in val.keys():
                val['organism'] = None

            # Convert to pdc_query
            print('INFO: Query dictionary passed to pseudomonas.com scraper will now be converted to a pdc_query object. See reference manual for more details.')
            val = _dict_to_pdc_query(**val)

        # Check keywords are internally consistent.
        if val.organism is not None and val.strain is not None:
            raise KeyError("Invalid combination of query keywords: 'organism' must not be combined with 'strain'.")

        # Check all values are strings or None.
        for v in val[:]:
            if not (isinstance(v, str) or v is None):
                raise TypeError("All values in the query must be of type str.")
        self.__query = val

    def connect(self):
        """ Establish a connection to pseudomonas.com by launching the selenium headless browser. """
        browser_opts = selenium.webdriver.firefox.options.Options()
        browser_opts.headless=True

        try:
            self.__browser = Browser(options=browser_opts)
            self.__browser.get(self.__pdc_url)

            print("INFO: Successfully connected to the pseudomonas.com database web frontend.")

        except:
            RuntimeError("Could not launch headless browser and connect to pseudomonas.com. Check your internet connection and mozilla firefox being installed.")

    def run_query(self):
        """ Submit a query to the db and get results as a structured table. Mutually exclusive with parameter 'organism'.
        """

        # Form http query string.
        _feature = self.query.feature
        if _feature is None:
            _feature = ''

        if self.query.strain is not None:
            _url = self.__pdc_url+"/primarySequenceFeature/list?c1=name&v1={0:s}&e1=1&term1={1:s}&assembly=complete".format(_feature, self.query.strain)
        elif query.organism is not None:
            _url = self.__pdc_url+"/primarySequenceFeature/list?c1=name&v1={0:s}&e1=1&term2={1:s}&assembly=complete".format(_feature, self.query.organism)

        print("DEBUG: Will now open {0:s}".format(_url))

        self.__browser.get(_url)

        # If we're looking for a unique feature.
        if _feature is not '':
            self.__browser.find_element_by_link_text(_feature.upper()).click()


        #tmp_tables = browser.find_element_by_xpath("/html/body/div[3]/div/div[3]/div[1]/")
        gene_overview_table = browser.find_element_by_xpath("/html/body/div[3]/div/div[3]/div[1]/table[1]")
        cross_references_table = browser.find_element_by_xpath("/html/body/div[3]/div/div[3]/div[1]/table[2]")
        product_table = browser.find_element_by_xpath("/html/body/div[3]/div/div[3]/div[2]/table[1]")
        subcellular_localization_table = browser.find_element_by_xpath("/html/body/div[3]/div/div[3]/div[2]/table[2]")
        pathogen_association_analysis_table = browser.find_element_by_xpath("/html/body/div[3]/div/div[4]/div[1]/table")
        comparative_genomics_table = browser.find_element_by_xpath("/html/body/div[3]/div/div[4]/div[2]/table")
        interactions_table = browser.find_element_by_xpath("/html/body/div[3]/div/div[5]/div[1]/table")
        references_table = browser.find_element_by_xpath("/html/body/div[3]/div/div[6]/div/div")

        self.__browser.save_screenshot(screenshot)

        self.__browser.quit()
        return screenshot

def _dict_to_pdc_query(**kwargs):
    """ """
    """
    Convert a dictionary of query key-value pairs to a pdc_query instance.

    :param kwargs: Dictionary of query key-value pairs.
    :type kwargs: dict
    """

    query = pdc_query(kwargs['strain'], kwargs['feature'], kwargs['organism'])

    return query

def pandasDF_from_heading(soup, table_heading):
    """ Find the table that belongs to the passed heading in a formatted html tree (the soup).

    :param soup: The html tree to parse.
    :type  soup: BeautifulSoup

    :param table_heading: The table heading to find.
    :type  table_heading: str

    :return: The table under the passed heading as a pandas.DataFrame
    :rtype: pandas.DataFrame

    """

    table_ht = str(soup.find(string=re.compile(table_heading)).find_next())
    #return table_ht

    return pandas.read_html(table_ht)[0]


