from setuptools import setup

setup(
    name='GenDBScraper',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    packages=['GenDBScraper',
              'GenDBScraper/Utilities'],
    license='MIT',
    long_description=open('README.md').read(),
)
