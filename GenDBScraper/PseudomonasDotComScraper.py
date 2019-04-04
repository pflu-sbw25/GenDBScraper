""" :module PseudomonasDotComScraper: Hosting the PseudomonasDotComScraper, an API for the https://www.pseudomonas.com database web interface. """

from GenDBScraper.Utilities.json_utilities import JSONEncoder
from GenDBScraper.Utilities.web_utilities import guarded_get, is_good_response

# 3rd party imports
from bs4 import BeautifulSoup
from collections import namedtuple
from doi2bib import crossref
from io import StringIO
from pubmed_lookup import Publication, PubMedLookup
import json
import logging
import numpy
import os
import pandas
import re
import tempfile

# Configure logging.
logging.basicConfig(format='%(asctime)s %(levelname)s: %(message)s', level=logging.INFO)

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
            val = [pdc_query(strain='sbw25')]


        exc = TypeError("The parameter 'query' must be a dict or pdc_query or a list, tuple, or set of queries. Examples: query={'strain' : 'sbw25', 'feature'='pflu0916'}; query=pdc_query(strain='sbw25', feature='pflu0916') or query=[pdc_query(strain='sbw25', feature='pflu0916'), pdc_query(strain='sbw25', feature='pflu0917')].")

        if not isinstance(val, list):
            if not (isinstance(val, dict) or isinstance(val, pdc_query)):
                raise exc
            else:
                val = [val]

        for i,v in enumerate(val):
            if isinstance(v, dict):
                pass
            elif isinstance(v, pdc_query):
                pass
            else:
                raise exc

        # Iterate over all queries.
        for i,v in enumerate(val):
            # Check keys if dict.
            if isinstance(v, dict):
                # Only these are acceptable query keywords.
                accepted_keys = ('strain', 'feature', 'organism')
                present_keys = v.keys()
                for k in present_keys:
                    if not k in accepted_keys:
                        raise KeyError("Only 'strain', 'feature', and 'organism' are acceptable keys.)")

                # Complete keywords.
                if not 'strain' in v.keys():
                    v['strain'] = None
                if not 'feature' in v.keys():
                    v['feature'] = None
                if not 'organism' in v.keys():
                    v['organism'] = None

                # Convert to pdc_query
                logging.info('Query dictionary passed to pseudomonas.com scraper will now be converted to a pdc_query object. See reference manual for more details.')
                v = _dict_to_pdc_query(**v)

            # Check keywords are internally consistent.
            if v.organism is not None and v.strain is not None:
                raise KeyError("Invalid combination of query keywords: 'organism' must not be combined with 'strain'.")

            # Check all values are strings or None.
            for vv in v[:]:
                if not (isinstance(vv, str) or vv is None):
                    raise TypeError("All values in the query must be of type str.")
            # Reset checked item.
            val[i] = v

        self.__query = val

    def connect(self):
        """ Connect to the database. """
        try:
            self.__browser = BeautifulSoup(guarded_get(self.__pdc_url), 'html.parser')
        except:
            self.__connected = False
            raise ConnectionError("Connecting to {0:s} failed. Make sure the URL is set correctly and is reachable.")

        self.__connected = True

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
        browser = BeautifulSoup(guarded_get(_url), 'html.parser')

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
        browser = BeautifulSoup(guarded_get(overview_url), 'lxml')

        # Loop over headings and get table as pandas.pandas.DataFrame.
        panels["Gene Feature Overview"] = _pandasDF_from_heading(browser, "Gene Feature Overview", None)

        # Get cross-references with
        self._get_cross_references(url, panels)

        #cr = _pandasDF_from_heading(browser, "Cross-References", None)
        #pandas.rename(index=int, columns={0: "type", 1 : "id"})
        #cr.iloc[:,"url"]
        #panels["Cross-References"] = cr

        panels["Product"] = _pandasDF_from_heading(browser, "Product", None)
        panels["Subcellular localization"] = _pandasDF_from_heading(browser, "Subcellular localization", None)
        panels["Pathogen Association Analysis"] = _pandasDF_from_heading(browser, "Pathogen Association Analysis", None)
        panels["Orthologs/Comparative Genomics"] = _pandasDF_from_heading(browser, "Orthologs/Comparative Genomics", None)
        panels["Interactions"] = _pandasDF_from_heading(browser, "Interactions", None)

        panels["References"] = _pandas_references(browser)

    def _get_cross_references(self, url, panels):
        """ Extract the cross-references table with hyperlinks from the feature overview tab. """
        # Get ovierview tab.
        operons_url = url + "&view=overview"
        soup = BeautifulSoup(guarded_get(operons_url), 'lxml')

        # Navigate to heading.
        table_heading = "Cross-References"
        heading = soup.find('h3', string=re.compile(table_heading))

        # Get content.
        cross_refs = heading.find_next_sibling('table')

        ref_types = []
        ref_ids = []
        ref_urls = []
        # Loop over operons.
        pattern = re.compile('[\t\s]')
        rows = cross_refs.find_all('tr')
        for i,row in enumerate(rows):
            cols = row.find_all('td')

            ref_type=cols[0].text
            ref_type = pattern.sub('', ref_type)

            hyperlink = cols[1].find_next('a')
            if hyperlink is not None:
                ref_id_text = pattern.sub('', hyperlink.text)
                ref_id_url = hyperlink.get('href')
            else:
                ref_id_text = cols[1].text
                ref_id_url = None

            ref_types.append(ref_type)
            ref_ids.append(ref_id_text)
            ref_urls.append(ref_id_url)

        df = pandas.DataFrame(numpy.empty((len(ref_types),0)))

        df.index = ref_types
        df['id'] = ref_ids
        df['url'] = ref_urls

        panels["Cross-References"] = df

    def _get_sequences(self, url, panels):
        """ Parse the 'Sequences' tab and extract the tables.

        :param url: The base URL of the feature.
        :type  url: str

        :param panels: The datastructure into which the tables are stored.
        :type  panel: dict

        """

        sequence_url = url +  "&view=sequence"
        browser = BeautifulSoup(guarded_get(sequence_url), 'html.parser')

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

        browser = BeautifulSoup(guarded_get(function_url), 'html.parser')

        panels["Gene Ontology"] = _pandasDF_from_heading(browser,"Gene Ontology", None)
        panels["Functional Classifications Manually Assigned by PseudoCAP"] = _pandasDF_from_heading(browser,"Functional Classifications Manually Assigned by PseudoCAP", None)
        panels["Functional Predictions from Interpro"] = _pandasDF_from_heading(browser,"Functional Predictions from Interpro", None)

        # Convert E-values to floats.
        panels["Functional Predictions from Interpro"]["E-value"] = pandas.to_numeric(panels["Functional Predictions from Interpro"]["E-value"], errors='coerce', downcast='float')

    def _get_motifs(self, url, panels):
        """ Parse the 'Motifs' tab and extract the tables.

        :param url: The base URL of the feature.
        :type  url: str

        :param panels: The datastructure into which the tables are stored.
        :type  panel: dict

        """

        # Get motifs tab.
        motifs_url = url + "&view=motifs"
        browser = BeautifulSoup(guarded_get(motifs_url), 'html.parser')

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
        soup = BeautifulSoup(guarded_get(operons_url), 'lxml')
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
        browser = BeautifulSoup(guarded_get(transposons_url), 'html.parser')

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
        browser = BeautifulSoup(guarded_get(updates_url), 'html.parser')

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
            doi_soup = BeautifulSoup(guarded_get(pubmed_link), 'lxml')
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


