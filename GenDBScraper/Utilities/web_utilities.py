""" :module: hosting various utilities built on top of the requests module. """

from contextlib import closing
import logging
from requests import get, post
from requests.exceptions import RequestException

def guarded_get(url):
    """ Get content of passed URL.

    :param url: The URL to parse.
    :type  url: str

    """
    # Safeguard opening the URL.
    with closing(get(url, stream=True, timeout=10)) as resp:
        if is_good_response(resp):
            logging.info("Connected to %s.", url)
            return resp.content
        else:
            raise RuntimeError("ERROR: Could not open "+url+" .")

def guarded_post(url, data):
    """ Post request to url in a safeguarded way. """

    try:
        resp = post(url, data=data, stream=True)
        if is_good_response(resp):
            logging.info("Connected to %s.", url)
            return resp
        else:
            raise RuntimeError("ERROR: Could not open "+url+" .")
    except:
        raise


def is_good_response(resp, expected_content_type='text'):
    """ Returns True if the response seems to be HTML, False otherwise.

    :param resp: The response to validate.
    :type  resp: http response as returned from contextlib.closing

    """

    content_type = resp.headers['Content-Type'].lower()
    return (resp.status_code == 200
            and content_type is not None
            and content_type.find(expected_content_type) > -1,
            )


