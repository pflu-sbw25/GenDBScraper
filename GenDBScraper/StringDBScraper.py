""" :module StringDBScraper: Hosting the StringDBScraper, an API for the https://string-db.org database web interface. """

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
        base_url = "http://string-db.org"
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

    def update_features(self):
        """ Replace the query features by the string-db identifiers. """
        resolved_ids = self.resolve_id(limit=1)
        self.query = stringdb_query(taxonId=self.query.taxonId, features=resolved_ids.preferredName.to_list())

    def resolve_id(self, **kwargs):
        """ Resolve the given identifier(s) to string-db.org's own identifiers.

        :param limit: (Optional): Limit the number of matches per query identifier (best matches come first). Default: limit=1
        :type  limit: int

        """
        """ Taken from  http://string-db.org/cgi/help.pl#Mapping-identifiers """

        if 'query' in kwargs.keys():
            self.query = kwargs['query']

        method = "get_string_ids"
        query_url = "/".join([self.base_url, 'api', 'json', method])

        data = dict(
                identifiers="\r".join(self.query.features),
                species    =self.query.taxonId if self.query.taxonId is not None else "",
                limit      =1 if not "limit" in kwargs.keys() else kwargs['limit'],
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

    def network_image(self, query=None, image_format='png', flavor=None, white_nodes=None, color_nodes=None, show_image=False):
        """ Grab the protein network image for given proteins (genes).

        :param query:  The (updated) query to submit.
        :type  query: (stringdb_query | dict)

        :param image_format: The image format for the network image (png, svg, hires_png)
        :type  image_format: str

        :param flavor: The type of network to draw between nodes (evidence, confidence (default), or actions).
        :type  flavor: str

        :param white_nodes: The number of white nodes to add. Default is 10 for single queries, 0 for multiple queries.
        :type  white_nodes: int

        :param color_nodes: The number of color nodes to add. Default is 0.
        :type  color_nodes: int

        :param show_image: Whether to render the image (default False). WARNING: untested feature.
        :type  show_image: bool

        """
        """ Inspired by  http://string-db.org/cgi/help.pl#Getting-STRING-network-image """

        if query is not None:
            self.query = query

        if not self.connected:
            raise IOError("Not connected to string-db.org.")


        format_map = {
                'png' : 'image',
                'image' : 'image',
                'hires_png': 'highres_image',
                'highres_image' : 'highres_image',
                'svg'           : 'svg',
                }

        method = "network"
        query_url = "/".join([self.base_url, 'api', format_map[image_format], method])

        data = dict(
                identifiers             = "\r".join(self.query.features),
                species                 = self.query.taxonId if self.query.taxonId is not None else "",
                add_white_nodes         = white_nodes,
                add_color_nodes         = color_nodes,
                required_score          = None,
                network_flavor          = None,
                caller_identity="https://gendbscraper.readthedocs.io",
                )

        # Get the response from post.
        response = web_utilities.guarded_post(query_url, data=data)


        # Determine file extension.
        suffix = ".png" if image_format.find("png") else "svg"

        # Setup image file.
        _, image_file = tempfile.mkstemp(prefix='string-db_network_', suffix=suffix)

        with open(image_file, 'wb') as image_fp:
            image_fp.write(response.content)

        if show_image:
            from PIL import Image
            Image.open(image_file).show()

        return image_file

    def network_interactions(self, nodes=None):
        """ Get the string-db network interactions as a pandas.DataFrame.

        :param nodes: The number of nodes to to add to the network based on their confidence score.
        :type  nodes: int

        """

        method = "network"
        query_url = "/".join([self.base_url, 'api', 'json', method])

        data = dict(
                identifiers             = "\r".join(self.query.features),
                species                 = self.query.taxonId if self.query.taxonId is not None else "",
                add_nodes               = nodes,
                required_score          = None,
                caller_identity="https://gendbscraper.readthedocs.io",
                )

        # Get the response from post.
        response = web_utilities.guarded_post(query_url, data=data)

        ret = pandas.DataFrame(response.json())

        return ret.reindex(columns = [
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
                )

    def interaction_partners(self, required_score=None, limit=None):
        """ Get the interaction partners.

        :param required_score: The minimum score for an interaction to be considered.
        :type  nodes: float

        :param limit: Limit the number of matches per query identifier (best matches come first). Default: limit=1
        :type  limit: int

        """

        method = "interaction_partners"
        query_url = "/".join([self.base_url, 'api', 'json', method])

        if limit is not None and not isinstance(limit, int):
            raise TypeError("limit must be an integer, {} was supplied.".format(type(limit)))
        if not isinstance(required_score, int):
            raise TypeError("required_score must be an integer (0 <= required_score <= 1000). It will be devided by 1000 to yield the actual minimum score cutoff.")

        data = dict(
                identifiers             = "\r".join(self.query.features),
                species                 = self.query.taxonId if self.query.taxonId is not None else "",
                required_score          = required_score,
                limit                   = limit,
                caller_identity="https://gendbscraper.readthedocs.io",
                )

        # Get the response from post.
        response = web_utilities.guarded_post(query_url, data=data)

        ret = pandas.DataFrame(response.json())

        return ret.reindex(
                columns = [
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
                )

    def similarity_scores(self):
        """ Get the interaction partners.

        :param required_score: The minimum score for an interaction to be considered.
        :type  nodes: float

        """

        raise NotImplementedError("This feature is currently not supported by string-db.org.")


        method = "homology"
        query_url = "/".join([self.base_url, 'api', 'json', method])

        data = dict(
                identifiers             = "\r".join(self.query.features),
                species                 = self.query.taxonId if self.query.taxonId is not None else "",
                caller_identity="https://gendbscraper.readthedocs.io",
                )

        # Get the response from post.
        response = web_utilities.guarded_post(query_url, data=data)

        ret = pandas.DataFrame(response.json())

        return ret.reindex(
                columns = [
                    'stringId_A',
                    'stringId_B',
                    'bitscore',
                    'start_A',
                    'end_A',
                    'start_B',
                    'end_B',
                    'size_B',
                    ]
                )

    def functional_enrichments(self):
        """ Get the interaction partners.

        """

        method = "enrichment"
        query_url = "/".join([self.base_url, 'api', 'json', method])

        data = dict(
                identifiers             = "\r".join(self.query.features),
                background_string_ids   = None,
                species                 = self.query.taxonId if self.query.taxonId is not None else "",
                caller_identity="https://gendbscraper.readthedocs.io",
                )

        # Get the response from post.
        response = web_utilities.guarded_post(query_url, data=data)

        # Setup and return dataframe.
        ret = pandas.DataFrame(response.json())

        return ret.reindex(
                columns = [
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
                )

    def interaction_enrichments(self, required_score=None):
        """ Get the interaction enrichments.

        :param required_score: The minimum score for an interaction to be considered.
        :type  nodes: float

        """

        method = "ppi_enrichment"
        query_url = "/".join([self.base_url, 'api', 'json', method])

        if required_score is not None and not isinstance(required_score, int):
            raise TypeError("required_score must be an integer (0 <= required_score <= 1000). It will be devided by 1000 to yield the actual minimum score cutoff.")

        data = dict(
                identifiers             = "\r".join(self.query.features),
                background_string_ids   = None,
                species                 = self.query.taxonId if self.query.taxonId is not None else "",
                caller_identity="https://gendbscraper.readthedocs.io",
                )

        # Get the response from post.
        response = web_utilities.guarded_post(query_url, data=data)

        # Setup and return dataframe.
        ret = pandas.DataFrame(response.json())

        return ret.reindex(
                columns = [
                    'number_of_nodes',
                    'number_of_edges',
                    'average_node_degree',
                    'local_clustering_coefficient',
                    'expected_number_of_edges',
                    'p_value',
                    ]
                )


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


