#!/usr/bin/env python

USERS={'admin': 'admin', 
        'user': 'user',
        'test': 'eW91IGFyZSBjcmF6eQo='}

"""
Twisted SSL webserver with basic authentication using plain in-memory passwords. 
The first argument is the path of the directory to serve; if not provided then the current folder is used (".").

INSTALL DEPENDENCIES:
    pip install twisted pyOpenSSL service_identity

GENERATE SSL CERTIFICATES:
    mkdir ~/.ssl && cd ~/.ssl
    openssl genrsa > privkey.pem
    openssl req -new -x509 -key privkey.pem -out cacert.pem -days 9999

USAGE:
    Requires running as root (normal users cannot bind to ports below 1024); 
    login with test_user/test_password

    sudo python twisted-web-ssl.py     # serve the current folder
    sudo python twisted-web-ssl.py /home
"""
import os
import sys

from twisted.web.static import File
from zope.interface import implements
from twisted.python import log
from twisted.internet import reactor, ssl
from twisted.web import server, resource, guard
from twisted.cred.portal import IRealm, Portal
from twisted.cred.checkers import InMemoryUsernamePasswordDatabaseDontUse
from twisted.python.log import startLogging

startLogging(sys.stdout)
home_dir = os.path.expanduser("~")

sslContext = ssl.DefaultOpenSSLContextFactory(
    'key.pem',
    'cert.pem'
)

class SimpleRealm(object):
    implements(IRealm)

    def __init__(self, path):
        self.path = path

    def requestAvatar(self, avatarId, mind, *interfaces):

        if resource.IResource in interfaces:
            return resource.IResource, File(self.path), lambda: None

        raise NotImplementedError()


def main(root):
    log.startLogging(sys.stdout)
    # alternative credential storage implementations: https://twistedmatrix.com/documents/current/api/twisted.cred.checkers.ICredentialsChecker.html
    checkers = [InMemoryUsernamePasswordDatabaseDontUse(**USERS)]

    wrapper = guard.HTTPAuthSessionWrapper(
        Portal(SimpleRealm(root), checkers),
        [guard.DigestCredentialFactory('md5', 'whatever.com')])

    reactor.listenSSL(443, server.Site(resource=wrapper),
                      contextFactory=sslContext)
    reactor.run()


if __name__ == '__main__':
    root = sys.argv[1] if len(sys.argv) > 1  else '.'
    main(root)