""" :module StringDBScraper: Hosting the StringDBScraper, an API for the https://www.pseudomonas.com database web interface. """

from GenDBScraper.RESTScraper import RESTScraper
from GenDBScraper.Utilities import web_utilities

# 3rd party imports
from collections import namedtuple
from doi2bib import crossref
from io import StringIO
from pubmed_lookup import Publication, PubMedLookup
import json
import logging
import os
import pandas
import re
import tempfile

# Configure logging.
logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO)

# Define the query datastructure.
stringdb_query = namedtuple(
        'stringdb_query',
        field_names=('taxonId', 'features'),
        defaults=('216595', []),
        )

class StringDBScraper(RESTScraper):
    """  An API for the string-db.org protein interaction database. """

    # Class constructor
    def __init__(self, query=None):
        """
        StringDBScraper constructor.

        :param query: The query to submit to string-db.org
        :type  query: (dict |
        """

        # Base class initialization.
        base_url = "https://string-db.org"
        super().__init__(base_url)

        self.query = query

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
        :type  val: (stringdb_query | dict)

        """

        # Checks
        if val is None:
            val = stringdb_query(taxonId=None, features=[])

        exc = TypeError("The parameter 'query' must be a dict or stringdb_query. Examples: query={'taxonId' : '216595', 'features'=['pflu0916']}; query=straindb_query(taxonId='216595', features['pflu0916', 'pflu0917']).")

        if not isinstance(val, (dict, stringdb_query)):
            raise exc

        # Check keys if dict.
        if isinstance(val, dict):
            # Only these are acceptable query keywords.
            accepted_keys = ('taxonId', 'features')
            present_keys = val.keys()
            for k in present_keys:
                if not k in accepted_keys:
                    raise KeyError("Only {0:s} are acceptable keys.".format(",".join(accepted_keys)))

            # Complete keywords.
            if not 'taxonId' in val.keys():
                val['taxonId'] = None
            if not isinstance(val['taxonId'], (str,int)):
                raise TypeError("taxonId must be a valid NCBI taxonId (str or int).")
            if not 'features' in val.keys():
                raise KeyError("You must specify a list of genes or products ('features').")

            # Convert to stringdb_query
            logging.info('Query dictionary passed to string-db scraper will now be converted to a stringdb_query object. See reference manual for more details.')

            val = stringdb_query(taxonId=val['taxonId'], features=val['features'])

        self.__query = val

    def resolve_id(self, **kwargs):
        """ Resolve the given identifier(s) to string-db.org's own identifiers.

        :param limit: (Optional): Limit the number of matches per query identifier (best matches come first). Default: limit=1
        :type  limit: int

        """
        """ Taken from  http://string-db.org/cgi/help.pl#Mapping-identifiers """

        method = "get_string_ids"
        query_url = "/".join([self.base_url, 'api', 'json', method])

        data = dict(
                identifiers="\r".join(self.query.features),
                species    =self.query.taxonId if self.query.taxonId is not None else "",
                limit      =1 if not "limit" in kwargs.keys() else limit,
                echo_query =1,
                caller_identity="https://gendbscraper.readthedocs.io",
                )

        # Get the response from post.
        response = web_utilities.guarded_post(query_url, data)

        ret = pandas.DataFrame(response.json())
        ret.index = ret['queryItem']
        del ret['queryItem']

        # Re-index.
        return ret.reindex(columns=['queryIndex', 'preferredName', 'stringId', 'ncbiTaxonId', 'taxonName', 'annotation'])

    def run_query(self, query=None):
        """ Run a query on pseudomonas.com

        :param query: The query object to run.
        :type  query: [list of] (pdc_query | dict)

        :return: The query results as a dictionary with 'strain_feature' keys.
        :rtype: dict

        """

        # Check if we're connected. Bail out if not.
        if not self.__connected:
            raise RuntimeError("Not connected. Call .connect() before submitting the query.")

        # If provided, update the local query object. This way, user can submit a query at run time.
        if query is not None:
            self.query = query

        results = dict()

        for query in self.query:
            key = "{0:s}__{1:s}".format(query.strain, query.feature)
            results[key] = self._run_one_query(query)

        return results

    def _get_feature_url(self, query):
        """ Get the base URL for the queried feature (gene).

        :param query: Query object.
        :type  query: pdc_query
        """

        # Form http self.query string.
        _feature = query.feature
        if _feature is None:
            _feature = ''

        # Assemble the html query.
        if query.strain is not None: # Searching for specific strain.
            _url = self.__pdc_url+"/primarySequenceFeature/list?c1=name&v1={0:s}&e1=1&term1={1:s}&assembly=complete".format(_feature, query.strain)
        elif query.organism is not None: # Searching for organism.
            _url = self.__pdc_url+"/primarySequenceFeature/list?c1=name&v1={0:s}&e1=1&term2={1:s}&assembly=complete".format(_feature, self.query.organism)

        # Debug info.
        logging.debug("Will now open {0:s} .".format(_url))

        # Get the soup for the assembled url.
        browser = BeautifulSoup(_guarded_get(_url), 'html.parser')

        # If we're looking for a unique feature.
        if _feature is not '':
            feature_link = browser.find_all('a', string=re.compile(_feature.upper()))[0].get('href')

        return self.__pdc_url + feature_link

    def _run_one_query(self, query):
        """ """
        """ Workhorse function to run a query.

        :param query: Query object to submit.
        :type  query: pdc_query
        """

        # Setup dict to store self.query results.
        panels = dict()
        feature_url =  self._get_feature_url(query)

        # Go through all panels and pull data.
        self._get_overview(feature_url, panels)
        self._get_sequences(feature_url, panels)
        self._get_functions_pathways_go(feature_url, panels)
        self._get_motifs(feature_url, panels)
        self._get_operons(feature_url, panels)
        self._get_transposon_insertions(feature_url, panels)
        self._get_updates(feature_url, panels)
        self._get_orthologs(feature_url, panels)

        # All done, return.
        return panels

    def _get_overview(self, url, panels):
        """ Parse the 'Overview' tab and extract the tables.

        :param url:  The base URL feature.
        :type  url: str

        :param panels [in/out]: The datastructure into which the tables are stored.
        :type  panel: dict

        """
        # Get overview data.
        overview_url = url + "&view=overview"

        # Get the soup.
        browser = BeautifulSoup(_guarded_get(overview_url), 'lxml')

        # Loop over headings and get table as pandas.pandas.DataFrame.
        panels["Gene Feature Overview"] = _pandasDF_from_heading(browser, "Gene Feature Overview", 0)
        panels["Cross-References"] = _pandasDF_from_heading(browser, "Cross-References", 0)
        panels["Product"] = _pandasDF_from_heading(browser, "Product", 0)
        panels["Subcellular localization"] = _pandasDF_from_heading(browser, "Subcellular localization",  0)
        panels["Pathogen Association Analysis"] = _pandasDF_from_heading(browser, "Pathogen Association Analysis", 0)
        panels["Orthologs/Comparative Genomics"] = _pandasDF_from_heading(browser, "Orthologs/Comparative Genomics", 0 )
        panels["Interactions"] = _pandasDF_from_heading(browser, "Interactions", 0)

        panels["References"] = _pandas_references(browser)

    def _get_sequences(self, url, panels):
        """ Parse the 'Sequences' tab and extract the tables.

        :param url: The base URL of the feature.
        :type  url: str

        :param panels: The datastructure into which the tables are stored.
        :type  panel: dict

        """

        sequence_url = url +  "&view=sequence"
        browser = BeautifulSoup(_guarded_get(sequence_url), 'html.parser')

        panels['Sequence Data'] = _pandasDF_from_heading(browser, "Sequence Data", None)

    def _get_functions_pathways_go(self, url, panels):
        """ Parse the 'Function/Pathways/GO' tab and extract the tables.

        :param url: The base URL of the feature.
        :type  url: str

        :param panels: The datastructure into which the tables are stored.
        :type  panel: dict

        """

        # Get functions, pathways, GO
        function_url = url + "&view=functions"

        browser = BeautifulSoup(_guarded_get(function_url), 'html.parser')

        panels["Gene Ontology"] = _pandasDF_from_heading(browser,"Gene Ontology", None)
        panels["Functional Classifications Manually Assigned by PseudoCAP"] = _pandasDF_from_heading(browser,"Functional Classifications Manually Assigned by PseudoCAP", None)
        panels["Functional Predictions from Interpro"] = _pandasDF_from_heading(browser,"Functional Predictions from Interpro", None)

    def _get_motifs(self, url, panels):
        """ Parse the 'Motifs' tab and extract the tables.

        :param url: The base URL of the feature.
        :type  url: str

        :param panels: The datastructure into which the tables are stored.
        :type  panel: dict

        """

        # Get motifs tab.
        motifs_url = url + "&view=motifs"
        browser = BeautifulSoup(_guarded_get(motifs_url), 'html.parser')

        panels["Motifs"] = None

    def _get_operons(self, url, panels):
        """ Parse the 'operons' tab and extract the tables.

        :param url: The base URL of the feature.
        :type  url: str

        :param panels: The datastructure into which the tables are stored.
        :type  panels: dict

        """

        # Get operons tab.
        operons_url = url + "&view=operons"
        soup = BeautifulSoup(_guarded_get(operons_url), 'lxml')
        table_heading = "Operons"

        # Navigate to heading.
        heading = soup.find('h3', string=re.compile(table_heading))

        # Get content.
        operons = heading.find_next_siblings('table')

        # Setup empty dict to store results.
        operons_dict = dict()

        # Loop over operons.
        for operon in operons:

            operon_dict = dict()

            try:
                tmp = pandas.read_html(str(operon))
            except:
                logging.warning("No operon data found.")
                break


            name = operon.findChild(string=re.compile("Operon name"))
            tabs = re.compile("\t*")
            name = tabs.sub("", name)
            name = name.split("\n")[2]

            operon_dict['Name'] = name
            operon_dict['Genes'] = tmp[1]

            evidence = str(operon.find(string=re.compile('Evidence')).find_next('div').text)
            evidence=re.compile("[\t\n\s\.]").sub("",evidence)
            evidence = re.sub("\.","",evidence)
            operon_dict['Evidence'] = evidence

            references = operon.find_all(string=re.compile('PubMed ID'))
            refs = []

            for ref in references:
                pubmed = ref.find_next_sibling('a')
                pubmed_url = pubmed.get('href')
                pubmed_id = str(pubmed.text)
                pubmed_id = re.compile('[\t\n\s]').sub('',pubmed_id)

                lookup = PubMedLookup(pubmed_id, '')
                citation = Publication(lookup).cite()

                refs.append(dict(pubmed_url=pubmed_url, citation=citation))
            operon_dict['References'] = pandas.DataFrame(refs)

            cross_references = str(operon.find(string=re.compile("Cross-References")).find_next('div').find_next('div').text)
            cross_references=re.compile("[\t\n\s]").sub("", cross_references)
            operon_dict['Cross-References'] = cross_references

            operons_dict[name] = operon_dict

        # Loop over headings and get table as pandas.pandas.DataFrame.
        panels['Operons'] = operons_dict

    def _get_transposon_insertions(self, url, panels):
        """ Parse the 'transposons' tab and extract the tables.

        :param url: The base URL of the feature.
        :type  url: str

        :param panels: The datastructure into which the tables are stored.
        :type  panel: dict

        """

        # Get transposons tab.
        transposons_url = url + "&view=transposons"
        browser = BeautifulSoup(_guarded_get(transposons_url), 'html.parser')

        table_heading = "Transposon Insertions"

        # Navigate to heading.
        headings = browser.find_all('h3', string=re.compile(table_heading))
        transposon_dict = dict()
        for h in headings:
            parent = h.parent
            try:
                td = pandas.read_html(str(parent))
            except ValueError:
                td = [pandas.DataFrame()]
                logging.warning("No table found, will return empty DataFrame.")
            key = h.get_text()
            key = re.compile("[\n\t]").sub("", key)

            transposon_dict[key] = td

        panels[table_heading] = transposon_dict

    def _get_updates(self, url, panels):
        """ Parse the 'Updates' tab and extract the tables.

        :param url: The base URL of the feature.
        :type  url: str

        :param panels: The datastructure into which the tables are stored.
        :type  panel: dict

        """

        # Get updates tab.
        updates_url = url + "&view=updates"
        browser = BeautifulSoup(_guarded_get(updates_url), 'html.parser')

        heading = browser.find('h3', string=re.compile('Annotation Updates'))
        annotation_table = pandas.read_html(str(heading.parent))

        panels['Annotation Updates'] = annotation_table

    def _get_orthologs(self, url, panels):
        """ Parse the 'Orthologs' tab and extract the tables.

        :param url: The base URL of the feature.
        :type  url: str

        :param panels: The datastructure into which the tables are stored.
        :type  panel: dict

        """

        ###
        # Orthologs are different in that they are queried from the pdc
        # orthologs database. We construct the corresponding URL and pull
        # the tab file directly.
        # TODO: or should we pull the fasta?

        # Get the pseudomonas.com id for this feature.
        pdc_id = url.split('id=')[1]

        # Construct the URL for the orthologs DB.
        orthologs_url = '/'.join([self.__pdc_url, 'orthologs', 'list?format=tab&extension=tab&id={}'.format(pdc_id)])

        # GET html. Bail out if none.
        try:
            request = get(orthologs_url)

            # Buffer the data.
            with StringIO(request.text) as stream:
                df = pandas.read_csv(stream, sep='\t')
                stream.close()

        except:
            logging.warning("No orthologs found. Will return empty DataFrame.")
            df = pandas.DataFrame()

        panels["Orthologs"] = df

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
            file_path = tempfile.mkstemp(prefix="pseudomonas_dot_com_query_", suffix=".json")[1]

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
    for query, result in loaded.items():
        ret[query] = {}
        for key, value in result.items():
            ret[query][key] = pandas.read_json(value)

    return ret

def _guarded_get(url):
    """ """
    """ Get content of passed URL to pass on to BeautifulSoup.

    :param url: The URL to parse.
    :type  url: str

    """

    # Safeguard opening the URL.
    with closing(get(url, stream=True, timeout=10)) as resp:
        if _is_good_response(resp):
            logging.info("Connected to %s.", url)
            return resp.content
        else:
            raise RuntimeError("ERROR: Could not open "+url+" .")

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
        logging.warning("No data found for %s. Will return empty pandas.DataFrame.", table_heading)
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
        # Get the link text.
        pubmed_link=a.get('href')

        citation = Publication(PubMedLookup(pubmed_link, '')).cite()
        raw.append(dict(pubmed_url=pubmed_link, citation=citation))

    # Return as pandas.DataFrame.
    return pandas.DataFrame(raw)

def _get_doi_from_ncbi(pubmed_link):
        """ Extract the DOI from a pubmed link. """

        if (pubmed_link != ''):
            doi_soup = BeautifulSoup(_guarded_get(pubmed_link), 'lxml')
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
    scraper = StringDBScraper(query)

    try:
        scraper.connect()
    except:
        logging.error("Could not connect to pseudomonas.com .")
        return 0

    # Run the query and serialize.
    try:
        results = scraper.run_query()
    except:
        logging.error("Query failed.")
        return 0

    try:
        path = scraper.to_json(results, args.outfile)
    except:
        logging.error("Could not write results to disk.")
        raise
        return 0

    # Message.
    logging.info("Query was successfull. Results stored in %s.", path)

    del scraper

    return 1

def _cleanup_str(chars, string):
    """ Replace all characters in chars in string by  "". """

    patterns = []
    for c in chars:
        patterns.append(re.compile(c))


    space = re.compile(" *")
    tab = re.compile("\t*")
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


