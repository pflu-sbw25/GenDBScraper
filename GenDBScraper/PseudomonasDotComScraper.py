""" :module PseudomonasDotComScraper: Hosting the PseudomonasDotComScraper, an API for the https://www.pseudomonas.com database web interface. """

from GenDBScraper.Utilities.json_utilities import JSONEncoder

# 3rd party imports
from bs4 import BeautifulSoup
from collections import namedtuple
from contextlib import closing
from doi2bib import crossref
from pubmed_lookup import Publication, PubMedLookup
from requests import get
from requests.exceptions import RequestException
import json
import os
import pandas
import re
import tempfile

# Define the query datastructure.
pdc_query = namedtuple('pdc_query',
        field_names=('strain', 'feature', 'organism'),
        defaults=(None, None, None),
        )

class PseudomonasDotComScraper():
    """  An API for the pseudomonas.com genome database using web scraping technology. """

    # Class constructor
    def __init__(self,
            query=None,
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
        self.__connected = False

        # Set attributes via setter.
        self.query = query

    # Attribute accessors
    @property
    def connected(self):
        return self.__connected

    @property
    def query(self):
        """ Get the query.

        :return: The query object.
        :rtype:  pdc_query

        """

        return self.__query

    @query.setter
    def query(self, val):
        """"""
        """ Set the query attribute.

        :param val: The value to set.
        :type  val: (pdc_query | dict)

        :raises KeyError: Both 'strain' and 'organism' are provided.
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
        """ Connect to the database. """
        try:
            self.__browser = BeautifulSoup(_simple_get(self.__pdc_url), 'html.parser')
        except:
            self.__connected = False
            raise ConnectionError("Connecting to {0:s} failed. Make sure the URL is set correctly and is reachable.")

        self.__connected = True

    def run_query(self, query=None):
        """ Submit a query to the db and get results.

        :param query: (Optional) the query object to submit.
        :type  query: pdc_query
        """

        # Check if we're connected. Bail out if not.
        if not self.__connected:
            raise RuntimeError("Not connected. Call .connect() before submitting the query.")

        # If provided, update the local query object. This way, user can submit a query at run time.
        if query is not None:
            self.query = query

        # Form http self.query string.
        _feature = self.query.feature
        if _feature is None:
            _feature = ''

        # Assemble the html query.
        if self.query.strain is not None: # Searching for specific strain.
            _url = self.__pdc_url+"/primarySequenceFeature/list?c1=name&v1={0:s}&e1=1&term1={1:s}&assembly=complete".format(_feature, self.query.strain)
        elif self.query.organism is not None: # Searching for organism.
            _url = self.__pdc_url+"/primarySequenceFeature/list?c1=name&v1={0:s}&e1=1&term2={1:s}&assembly=complete".format(_feature, self.query.organism)

        # Debug info.
        print("DEBUG: Will now open {0:s}".format(_url))

        # Get the soup for the assembled url.
        browser = BeautifulSoup(_simple_get(_url), 'html.parser')

        # If we're looking for a unique feature.
        if _feature is not '':
            feature_link = browser.find_all('a', string=re.compile(_feature.upper()))[0].get('href')

        # Prepend base url.
        feature_link = self.__pdc_url+feature_link
        # Get the soup.
        #browser = BeautifulSoup(_simple_get(feature_link), 'html.parser')
        browser = BeautifulSoup(_simple_get(feature_link), 'lxml')

        # Setup dict to store self.query results.
        panels = dict()

        # Loop over headings and get table as pandas.pandas.DataFrame.
        panels["Gene Feature Overview"] = _pandasDF_from_heading(browser, "Gene Feature Overview", 0)
        panels["Cross-References"] = _pandasDF_from_heading(browser, "Cross-References", 0)
        panels["Product"] = _pandasDF_from_heading(browser, "Product", 0)
        panels["Subcellular localization"] = _pandasDF_from_heading(browser, "Subcellular localization",  0)
        panels["Pathogen Association Analysis"] = _pandasDF_from_heading(browser, "Pathogen Association Analysis", 0)
        panels["Orthologs/Comparative Genomics"] = _pandasDF_from_heading(browser, "Orthologs/Comparative Genomics", 0 )
        panels["Interactions"] = _pandasDF_from_heading(browser, "Interactions", 0)

        panels["References"] = _pandas_references(browser)

        # Assemble url for functions (tab "Function/Pathways/GO")
        function_url = feature_link + "&view=functions"
        browser = BeautifulSoup(_simple_get(function_url), 'html.parser')

        panels["Gene Ontology"] = _pandasDF_from_heading(browser,"Gene Ontology", None)
        panels["Functional Classifications Manually Assigned by PseudoCAP"] = _pandasDF_from_heading(browser,"Functional Classifications Manually Assigned by PseudoCAP", None)
        panels["Functional Predictions from Interpro"] = _pandasDF_from_heading(browser,"Functional Predictions from Interpro", None)

        # Return.
        return panels

    def to_json(self, results, outfile=None):
        """ Serialize results dictionary to json.

        :param results: The results dictionary (dict of pandas.DataFrame).
        :type  results: dict

        :param outfile: Path to file for writing query results to. Default: None, will write to temp file.
        :type  outfile: str

        :raises IOError: 'outfile' not writable.

        :return: If successful, path to written file.
        """

        if outfile is None:
            # Setup the filename from query items.
            if self.query.strain is not None:
                major = self.query.strain
            else:
                major = self.query.organism
            minor = self.query.feature
            file_path = tempfile.mkstemp(prefix="{0:s}_{1:s}_".format(major, minor), suffix=".json")[1]

        else:
            file_path = outfile

        # Call the workhorse.
        _serialize(file_path, results)

        return file_path

    def from_json(self, infile):
        """ Deserialize a json file into a results dictionary (a dict of pandas.DataFrame).

        :param infile: The file path of the json file to load.
        :type  infile: str

        """

        return _deserialize(infile)

def _serialize(path, obj):
    """ """
    """ Serialize the passed dictionary (obj) to path. """

    with open(path, 'w') as fp:
        json.dump(obj, fp, cls=JSONEncoder)

def _deserialize(path):
    """ """
    """ Deserialize a json file (located at 'path') into a dictionary. Reconstruct pandas.DataFrames from loaded content. """
    with open(path, 'r') as fp:
        loaded = json.load(fp)

    ret = {}
    for k,v in loaded.items():
        ret[k] = pandas.read_json(v)

    return ret

def _simple_get(url):
    """ """
    """ Get content of passed URL to pass on to BeautifulSoup.

    :param url: The URL to parse.
    :type  url: str

    """

    # Safeguard opening the URL.
    with closing(get(url, stream=True)) as resp:
        if _is_good_response(resp):
            print("INFO: Good response from "+url+".")
            return resp.content
        else:
            raise RuntimeError("ERROR: Could not open "+url+".")

def _is_good_response(resp):
    """ """
    """ Returns True if the response seems to be HTML, False otherwise.

    :param resp: The response to validate.
    :type  resp: http response as returned from contextlib.closing

    """

    content_type = resp.headers['Content-Type'].lower()
    return (resp.status_code == 200
            and content_type is not None
            and content_type.find('html') > -1)

def _dict_to_pdc_query(**kwargs):
    """ """
    """
    Convert a dictionary of query key-value pairs to a pdc_query instance.

    :param kwargs: Dictionary of query key-value pairs.
    :type kwargs: dict
    """

    query = pdc_query(kwargs['strain'], kwargs['feature'], kwargs['organism'])

    return query

def _pandasDF_from_heading(soup, table_heading, index_column=0):
    """ """
    """ Find the table that belongs to the passed heading in a formatted html tree (the soup).

    :param soup: The html tree to parse.
    :type  soup: BeautifulSoup

    :param table_heading: The table heading to find.
    :type  table_heading: str

    :param index_column: Which column to use as the pandas.DataFrame's index.
    :type  index_column: int

    :return: The table under the passed heading as a pandas.pandas.DataFrame
    :rtype: pandas.pandas.DataFrame

    """

    # Get table html string.
    table_ht = str(soup.find('h3', string=re.compile(table_heading)).find_next())

    try:
        df = pandas.read_html(table_ht, index_col=index_column)[0]
    except:
        print("WARNING: No data found for '"+table_heading+"'. Will return empty pandas.DataFrame.")

        df = pandas.DataFrame()

    return df

def _pandas_references(soup):
    """ Extract references from given html soup and return them as pandas pandas.DataFrame. """

    # Setup container to store parsed information.
    raw = []

    # Get the References "table".
    ref_soup = soup.find("h3", string=re.compile('^References'))

    # Get all <a> tags.
    a_tags = ref_soup.find_next().find_all('a')

    # Loop over all <a> tags
    for i,a in enumerate(a_tags):
        #pubmed_id_string = a.string

        # Strip \t, \n, and spaces.
        #pubmed_id = re.sub(pattern="[\t,\n,\s]", repl="", string=pubmed_id_string)
        # Get the link text.
        pubmed_link=a.get('href')

        citation = Publication(PubMedLookup(pubmed_link, '')).cite()

        ## Append to storage container.


        ## Get doi from pubmed link.
        #doi =_get_doi_from_ncbi(pubmed_link)

        ## Get bibliographic information from doi.
        #bib = _get_bib_from_doi(doi)

        raw.append(dict(pubmed_url=pubmed_link, citation=citation))


    # Return as pandas.DataFrame.
    df = pandas.DataFrame(raw)

    return df

def _get_doi_from_ncbi(pubmed_link):
        """ Extract the DOI from a pubmed link. """

        if (pubmed_link != ''):
            doi_soup = BeautifulSoup(_simple_get(pubmed_link), 'lxml')
        line = doi_soup.find(string=re.compile("DOI")).find_parent().find_parent()
        a = line.find('a', string=re.compile('10\.[0-9]*\/'))
        doi_string = a.text
        doi = re.sub("[\t,\n,\s]","",doi_string)

        return doi

def _get_bib_from_doi(doi):
    """ Get bibliographic information from a given doi."""

    # Get bib data.
    success, json = crossref.get_json(doi)

    if success and json['status'].lower() == 'ok':
        message = json['message']

        entry = {'doi'      :doi,
                 'first_author'   :"{0:s}, {1:s}".format(message['author'][0]['family'],
                                                   message['author'][0]['given']),
                 'title'    :message['title'][0],
                 'container':message['container-title'][0],
                 'volume'   :message['volume'],
                 'page'     :message['page'],
                 'date'     :"{0:d}-{1:02d}-{2:02d}".format(*(message['published-print']['date-parts'][0])),
                }
    return entry

def _run_from_cli(args):
    """ Called if run via command line interface.

    :param args: Command line arguments.
    :type  args: argparse.ArgumentsObject

    """

    # Construct the query.
    query = pdc_query(args.strain, args.feature, args.organism)

    # Construct the Scraper.
    scraper = PseudomonasDotComScraper(query)

    try:
        scraper.connect()
    except:
        print("ERROR: Could not connect to pseudomonas.com.")
        return 0

    # Run the query and serialize.
    try:
        results = scraper.run_query()
    except:
        print("ERROR: Query failed.")
        return 0

    try:
        path = scraper.to_json(results, args.outfile)
    except:
        print("ERROR: Could not write results to disk.")
        raise
        return 0

    # Message.
    print("INFO: Query was successfull. Results stored in {0:s}.".format(path))

    del scraper

    return 1

if __name__ == "__main__":

    from argparse import ArgumentParser

    # Setup argument parser.
    parser = ArgumentParser()

    parser.add_argument("-o",
                        "--outfile",
                        dest="outfile",
                        default=None,
                        required=False,
                        help="Where to write the query results.",
                        )

    parser.add_argument("-f",
                        "--feature",
                        dest="feature",
                        default=None,
                        required=True,
                        help="The gene/feature to query from pseudomonas.com.")

    org_group = parser.add_mutually_exclusive_group(required=True)
    org_group.add_argument("-s",
                        "--strain",
                        dest="strain",
                        default=None,
                        help="The strain to query from pseudomonas.com. Mutually exclusive with parameter -o/--organism option.")

    org_group.add_argument("-O",
                        "--organism",
                        dest="organism",
                        default=None,
                        help="The organism to query from pseudomonas.com. Mutually exclusive with parameter 'strain'.")

    # Parse arguments.
    args = parser.parse_args()

    _run_from_cli(args)


